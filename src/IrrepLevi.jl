# ═══════════════════════════════════════════════════════════════════════════════
#  IrrepLevi — irreducible representation of the Levi subgroup
#
#  An irreducible representation of the Levi subgroup L of a parabolic P
#  is encoded by a P-dominant weight λ in the weight lattice of G.  It
#  can also be described as a character of the center Z(L) (the "central
#  part") tensored with an irreducible representation of the semisimple
#  part [L,L] (the "semisimple part").
#
#  The decomposition uses the decomposition_matrix change-of-basis from
#  MarkedDynkinType.jl.
# ═══════════════════════════════════════════════════════════════════════════════

export IrrepLevi
export central_part, semisimple_part
export to_ambient_weight, fiber_dimension, p_dominant_weight
export central_scaling_factor

# Names from Lie and StaticArrays are available via the parent module's
# `using Lie` and `using StaticArrays`.

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _central_length(::Type{MDT}) -> Int

Return the number of marked (non-parabolic) nodes, i.e. the compile-time
length `K` of the `central::SVector{K,Int}` field in `IrrepLevi{MDT,K}`.
"""
@generated function _central_length(
  ::Type{MDT},
) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  K = length(Marked)
  return :($K)
end

"""
    IrrepLevi{MDT}

An irreducible representation of the Levi subgroup associated to the
marked Dynkin type `MDT`, encoded by its P-dominant weight `λ` in the
weight lattice of the ambient group ``G``.

Internally the weight also stores the decomposition into:
- `central::SVector{K,Int}`: scaled coordinates of the central character,
  indexed by the marked (nonparabolic) nodes. Stored as integers
  multiplied by [`central_scaling_factor`](@ref) for efficiency.
- `semisimple::WeightLatticeElem`: highest weight of the semisimple part,
  as a weight in the Levi's weight lattice

The decomposition is via the [`decomposition_matrix`](@ref) change of basis.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{3}, 1);

julia> rep = IrrepLevi(MDT, ω₁);

julia> central_part(rep)
1-element Vector{Rational{Int64}}:
 1//2

julia> semisimple_part(rep)
ω1
```
"""
struct IrrepLevi{MDT<:MarkedDynkinType, K}
  λ::WeightLatticeElem             # the P-dominant weight (ambient weight lattice)
  central::SVector{K,Int}          # scaled central character (×central_scaling_factor)
  semisimple::WeightLatticeElem    # highest weight of semisimple part (Levi lattice)
end

# Convenience outer constructor: infers K from the SVector.
IrrepLevi{MDT}(λ, c::SVector{K,Int}, ss) where {MDT,K} = IrrepLevi{MDT,K}(λ, c, ss)

# ─── Core accessors ──────────────────────────────────────────────────────────

"""
    p_dominant_weight(rep::IrrepLevi{MDT}) -> WeightLatticeElem{DT,R}

Return the P-dominant weight ``\\lambda`` that defines this Levi representation.

This accessor is `@generated` to ensure the return type
`WeightLatticeElem{DT,R}` is fully inferred from `MDT`.
"""
@generated function p_dominant_weight(rep::IrrepLevi{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  R = rank(DT)
  return :(rep.λ::WeightLatticeElem{$DT,$R})
end

# ─── Central scaling factor ──────────────────────────────────────────────────

"""
    central_scaling_factor(::Type{MDT}) -> Int

Return the scaling factor for the central part of Levi representations
associated to `MDT`.  This is the LCM of all denominators of the
inverse Cartan matrix entries at the marked rows:

``\\mathrm{lcm}\\{\\mathrm{denom}(C^{-1}[j, k]) : j \\in \\mathrm{Marked},\\; 1 \\le k \\le R\\}``

Central characters are stored internally as integers multiplied by this
factor, eliminating all `Rational{Int}` arithmetic from hot paths.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> central_scaling_factor(MarkedDynkinType{TypeA{3}, (2,)})
4

julia> central_scaling_factor(MarkedDynkinType{TypeA{4}, (2,)})
5
```
"""
@generated function central_scaling_factor(
  ::Type{MDT},
) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  R = rank(DT)
  Cinv = Lie.cartan_matrix_inverse(DT)
  sf = 1
  for j in Marked
    for k in 1:R
      sf = lcm(sf, denominator(Cinv[j, k]))
    end
  end
  return :($sf)
end

"""
    _apply_central_ext(::Type{MDT}, λ_ivec::SVector{R,Int}) -> SVector{K,Int}

Apply the central-extraction map inline: emits
`result[i] = Σ_k (sf * C⁻¹[Marked[i],k]) * λ_ivec[k]`
as direct scalar arithmetic baked into the generated code,
so no `SMatrix` or `Rational` type appears at runtime.
"""
@generated function _apply_central_ext(
  ::Type{MDT},
  λ_ivec::SVector{R0,Int},
) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked,R0}
  R = rank(DT)
  Cinv = Lie.cartan_matrix_inverse(DT)
  K = length(Marked)
  sf = 1
  for j in Marked, k in 1:R
    sf = lcm(sf, denominator(Cinv[j, k]))
  end
  row_exprs = map(1:K) do i
    terms = [:($(round(Int, Cinv[Marked[i], k] * sf)) * λ_ivec[$k]) for k in 1:R]
    length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)
  end
  return :(SVector{$K,Int}($(row_exprs...)))
end

"""
    _amb_scalars(::Type{MDT}) -> (sf_total, ratio)

At code-generation time, compute the two integer scaling factors needed by
`IrrepLevi(MDT, central, semisimple)` to reconstruct the ambient weight:
- `sf_total`: lcm of `sf_central` and all denominators of `Minv`
- `ratio = sf_total ÷ sf_central`

Returns only plain integers (no `SMatrix`), so no StaticArrays specialization
is triggered by the return type.
"""
@generated function _amb_scalars(
  ::Type{MDT},
) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  R = rank(DT)
  Cinv = Lie.cartan_matrix_inverse(DT)
  sf = 1
  for j in Marked, k in 1:R
    sf = lcm(sf, denominator(Cinv[j, k]))
  end
  unmarked = [i for i in 1:R if !(i in Marked)]
  M = zeros(Rational{Int}, R, R)
  for i in unmarked
    M[i, i] = 1
  end
  for j in Marked, k in 1:R
    M[j, k] = Cinv[j, k]
  end
  Minv = inv(M)
  sf_total = sf
  for j in 1:R, k in 1:R
    sf_total = lcm(sf_total, denominator(Minv[j, k]))
  end
  ratio = sf_total ÷ sf
  return :(($(sf_total), $(ratio)))
end

"""
    _apply_Minv_int(::Type{MDT}, x::SVector{R,Int}) -> SVector{R,Int}

Apply `sf_total * Minv` inline: emits
`result[i] = Σ_j (sf_total * Minv[i,j]) * x[j]`
as direct scalar arithmetic baked into the generated code.
No `SMatrix` type is constructed or dispatched at runtime, avoiding
the StaticArrays `gen_by_access` / `mul_parent` specialization cost.
"""
@generated function _apply_Minv_int(
  ::Type{MDT},
  x::SVector{R0,Int},
) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked,R0}
  R = rank(DT)
  Cinv = Lie.cartan_matrix_inverse(DT)
  sf = 1
  for j in Marked, k in 1:R
    sf = lcm(sf, denominator(Cinv[j, k]))
  end
  unmarked = [i for i in 1:R if !(i in Marked)]
  M = zeros(Rational{Int}, R, R)
  for i in unmarked
    M[i, i] = 1
  end
  for j in Marked, k in 1:R
    M[j, k] = Cinv[j, k]
  end
  Minv = inv(M)
  sf_total = sf
  for j in 1:R, k in 1:R
    sf_total = lcm(sf_total, denominator(Minv[j, k]))
  end
  row_exprs = map(1:R) do i
    terms = [:($(round(Int, Minv[i, j] * sf_total)) * x[$j]) for j in 1:R]
    length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)
  end
  return :(SVector{$R,Int}($(row_exprs...)))
end

"""
    central_part(rep::IrrepLevi) -> Vector{Rational{Int}}

Return the central character part of the Levi representation.

Internally the central character is stored as scaled integers
(multiplied by [`central_scaling_factor`](@ref)).  This accessor
unscales and returns the original `Rational{Int}` values.
"""
function central_part(rep::IrrepLevi{MDT}) where {MDT}
  sf = central_scaling_factor(MDT)
  return Vector{Rational{Int}}(Rational{Int}[c // sf for c in rep.central])
end

"""
    semisimple_part(rep::IrrepLevi{MDT}) -> WeightLatticeElem{LT,LR}

Return the highest weight of the semisimple part.

This accessor is `@generated` to ensure the return type
`WeightLatticeElem{LT,LR}` is fully inferred from `MDT`,
eliminating dynamic dispatch on downstream Lie.jl operations
(`degree`, `dual`, `tensor_product`, `exterior_power`, etc.).
"""
@generated function semisimple_part(rep::IrrepLevi{MDT}) where {MDT<:MarkedDynkinType}
  LT = levi_type(MDT)
  if LT === nothing
    return :(rep.semisimple::WeightLatticeElem{TypeA{1},1})
  else
    LR = rank(LT)
    return :(rep.semisimple::WeightLatticeElem{$LT,$LR})
  end
end

# ─── Construction from ambient weight ───────────────────────────────────────

"""
    IrrepLevi(::Type{MDT}, λ::WeightLatticeElem) -> IrrepLevi{MDT}

Construct an irreducible Levi representation from a weight ``\\lambda``
of the ambient group ``G``, by applying the [`decomposition_matrix`](@ref)
change of basis to decompose into central + semisimple parts.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{4}, 1);

julia> rep = IrrepLevi(MDT, ω₁);

julia> fiber_dimension(rep)
2
```
"""
function IrrepLevi(::Type{MDT}, λ::WeightLatticeElem) where {
  MDT<:MarkedDynkinType
}
  LT = levi_type(MDT)
  unmarked = unmarked_nodes(MDT)

  # Use integer arithmetic throughout: the unmarked rows of the decomposition
  # matrix are identity rows, so new_coords[u] == λ[u] for every unmarked u.
  # The marked rows are sf-scaled integer multiples of C⁻¹; _apply_central_ext
  # bakes them in at compile time and emits the multiply as inline scalar
  # arithmetic — no SMatrix type appears at runtime.
  R = rank(_ambient_type(MDT))
  λ_ivec = SVector{R,Int}(Tuple(coefficients(λ)))
  central = _apply_central_ext(MDT, λ_ivec)  # SVector{K,Int} — no conversion needed

  # Semisimple part: read off unmarked coordinates directly from λ.
  # Build as SVector{LR,Int} to hit the SVector overload of WeightLatticeElem
  # and avoid constructing a heap-allocated Array.
  if LT === nothing
    semisimple = WeightLatticeElem(TypeA{1}, SVector{1,Int}(0))
  else
    LR = rank(LT)
    perm = levi_permutation(MDT)
    ss_coords = SVector{LR,Int}(ntuple(j -> λ_ivec[unmarked[perm[j]]], LR))
    semisimple = WeightLatticeElem(LT, ss_coords)
  end

  return IrrepLevi{MDT}(λ, central, semisimple)
end

"""
    IrrepLevi(::Type{MDT}, central::AbstractVector{Int}, semisimple) -> IrrepLevi{MDT}

Construct an `IrrepLevi` from its **scaled** central character `central`
(pre-multiplied by [`central_scaling_factor`](@ref)) and the highest weight
`semisimple` of its semisimple part.  Any `AbstractVector{Int}` is accepted
and converted to the internal `SVector{K,Int}` representation.
The ambient P-dominant weight ``\\lambda`` is recovered automatically.
"""
function IrrepLevi(::Type{MDT}, central::AbstractVector{Int}, semisimple::WeightLatticeElem) where {
  MDT<:MarkedDynkinType
}
  K = _central_length(MDT)
  IrrepLevi(MDT, SVector{K,Int}(central), semisimple)
end

function IrrepLevi(::Type{MDT}, central::SVector{K,Int}, semisimple::WeightLatticeElem) where {
  MDT<:MarkedDynkinType,K
}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  unmarked = unmarked_nodes(MDT)
  LT = levi_type(MDT)
  R = rank(DT)

  sf_total, ratio = _amb_scalars(MDT)

  # Build integer coordinate vector: scaled by sf_total.
  # Use plain Vector{Int} — no SizedArray wrapper, no abstract-iteration inference.
  #   coords_full[m] = central[idx] * ratio  (= sf_total/sf_central * central)
  #   coords_full[u] = ss_nat[u]    * sf_total
  coords_full = Vector{Int}(undef, R)
  for (idx, m) in enumerate(Marked)
    coords_full[m] = central[idx] * ratio
  end

  if length(unmarked) > 0
    ss_vec = Lie.coefficients(semisimple)
    LR = length(ss_vec)
    if LR > 0
      perm = levi_permutation(MDT)
      inv_perm = Vector{Int}(undef, LR)
      for j in 1:LR
        inv_perm[perm[j]] = j
      end
      for (i, u) in enumerate(unmarked)
        coords_full[u] = ss_vec[inv_perm[i]] * sf_total
      end
    else
      for u in unmarked
        coords_full[u] = 0
      end
    end
  end

  # _apply_Minv_int emits sf_total*Minv inline; result = sf_total² * λ.
  # map over SVector returns SVector, so WeightLatticeElem hits its SVector overload.
  ambient_scaled = _apply_Minv_int(MDT, SVector{R,Int}(coords_full))
  sf_sq = sf_total * sf_total
  λ = WeightLatticeElem(DT, map(c -> div(c, sf_sq), ambient_scaled))

  return IrrepLevi{MDT}(λ, central, semisimple)
end

"""
    IrrepLevi(::Type{MDT}, central::Vector{Rational{Int}}, semisimple) -> IrrepLevi{MDT}

Construct an `IrrepLevi` from its central character `central` (a
`Vector{Rational{Int}}` indexed by marked nodes) and the highest weight
`semisimple` of its semisimple part.  The rational values are converted
to scaled integers internally.
"""
function IrrepLevi(::Type{MDT}, central::Vector{Rational{Int}}, semisimple::WeightLatticeElem) where {
  MDT<:MarkedDynkinType
}
  sf = central_scaling_factor(MDT)
  central_scaled = Int[Int(c * sf) for c in central]
  IrrepLevi(MDT, central_scaled, semisimple)
end

# ─── Back-conversion to ambient weight ───────────────────────────────────────

"""
    to_ambient_weight(::Type{MDT}, rep::IrrepLevi{MDT}) -> WeightLatticeElem

Convert an `IrrepLevi` back to a weight in the ambient group's weight lattice.

Since the `IrrepLevi` already stores the ambient P-dominant weight ``\\lambda``,
this simply returns it via `p_dominant_weight`.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{3}, 1);

julia> rep = IrrepLevi(MDT, ω₁);

julia> to_ambient_weight(MDT, rep) == ω₁
true
```
"""
function to_ambient_weight(::Type{MDT}, rep::IrrepLevi{MDT}) where {
  MDT<:MarkedDynkinType
}
  p_dominant_weight(rep)
end

# ─── Fiber dimension ─────────────────────────────────────────────────────────

"""
    fiber_dimension(rep::IrrepLevi{MDT}) -> BigInt

Return the dimension of the fiber (= dimension of the irreducible
representation of the semisimple part of the Levi).

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{3}, 1);

julia> fiber_dimension(IrrepLevi(MDT, ω₁))
2

julia> ω₂ = fundamental_weight(TypeA{3}, 2);

julia> fiber_dimension(IrrepLevi(MDT, ω₂))
1
```
"""
function fiber_dimension(rep::IrrepLevi{MDT}) where {MDT}
  LT = levi_type(MDT)
  LT === nothing && return BigInt(1)
  ss = semisimple_part(rep)
  iszero(ss) && return BigInt(1)
  is_dominant(ss) || return BigInt(0)
  degree(ss)
end

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, rep::IrrepLevi{MDT}) where {MDT}
  # Print the P-dominant weight directly
  print(io, "(", sprint(show, p_dominant_weight(rep)), ")")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Monoidal operations
# ═══════════════════════════════════════════════════════════════════════════════

"""
    tensor_product(a::IrrepLevi{MDT}, b::IrrepLevi{MDT}) -> Vector{IrrepLevi{MDT}}

Compute the tensor product of two irreducible Levi representations.

The central parts add, and the semisimple parts tensor (decomposing into
irreducibles via the Weyl character ring).

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{3}, 1);

julia> ω₃ = fundamental_weight(TypeA{3}, 3);

julia> reps_a = IrrepLevi(MDT, ω₁);

julia> reps_b = IrrepLevi(MDT, ω₃);

julia> result = tensor_product(reps_a, reps_b);

julia> length(result)
1
```
"""
function tensor_product(a::IrrepLevi{MDT,K}, b::IrrepLevi{MDT,K}) where {MDT,K}
  LT = levi_type(MDT)

  # Central parts add (in scaled representation); SVector addition is type-stable.
  new_central = a.central + b.central

  # Semisimple parts: tensor product decomposition
  ss_a = semisimple_part(a)
  ss_b = semisimple_part(b)
  if LT === nothing || (iszero(ss_a) && iszero(ss_b))
    return [IrrepLevi(MDT, new_central, ss_a)]
  end

  if iszero(ss_a)
    return [IrrepLevi(MDT, new_central, ss_b)]
  end
  if iszero(ss_b)
    return [IrrepLevi(MDT, new_central, ss_a)]
  end

  # Tensor product via Lie.jl's Weyl character ring
  χ = Lie.tensor_product(ss_a, ss_b)

  result = Vector{IrrepLevi{MDT,K}}()
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(MDT, new_central, hw))
    end
  end
  return result
end

"""
    dual(rep::IrrepLevi{MDT}) -> IrrepLevi{MDT}

Compute the dual of an irreducible Levi representation.
The central part is negated, and the semisimple part is dualized
(via the action of the longest Weyl group element).

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{3}, 1);

julia> rep = IrrepLevi(MDT, ω₁);

julia> d = dual(rep);

julia> central_part(d) == -central_part(rep)
true
```
"""
function dual(rep::IrrepLevi{MDT}) where {MDT}
  LT = levi_type(MDT)
  new_central = -rep.central
  ss = semisimple_part(rep)

  if LT === nothing || iszero(ss)
    return IrrepLevi(MDT, new_central, ss)
  end

  new_ss = Lie.dual(ss)
  return IrrepLevi(MDT, new_central, new_ss)
end

"""
    exterior_power(rep::IrrepLevi{MDT}, k::Int) -> Vector{IrrepLevi{MDT}}

Compute the k-th exterior power of an irreducible Levi representation.
The central part scales by k, and the semisimple part uses ⋀ᵏ.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (3,)};

julia> ω₁ = fundamental_weight(TypeA{4}, 1);

julia> rep = IrrepLevi(MDT, ω₁);

julia> result = exterior_power(rep, 0);

julia> length(result)
1
```
"""
function exterior_power(rep::IrrepLevi{MDT,K}, k::Integer) where {MDT,K}
  k = Int(k)
  LT = levi_type(MDT)

  k < 0 && return Vector{IrrepLevi{MDT,K}}()
  k == 0 && return [IrrepLevi(MDT, zero(rep.central),
                              WeightLatticeElem(LT === nothing ? TypeA{1} : LT,
                                               zeros(Int, rank(LT === nothing ? TypeA{1} : LT))))]
  k == 1 && return [rep]

  new_central = k * rep.central
  ss = semisimple_part(rep)

  if LT === nothing || iszero(ss)
    return Vector{IrrepLevi{MDT,K}}()  # ⋀^k of a 1-dim rep for k > 1 is 0
  end

  dim_ss = Int(degree(ss))
  k > dim_ss && return Vector{IrrepLevi{MDT,K}}()

  χ = Lie.exterior_power(ss, k)

  result = Vector{IrrepLevi{MDT,K}}()
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(MDT, new_central, hw))
    end
  end
  return result
end

"""
    symmetric_power(rep::IrrepLevi{MDT}, k::Int) -> Vector{IrrepLevi{MDT}}

Compute the k-th symmetric power of an irreducible Levi representation.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{3}, 1);

julia> rep = IrrepLevi(MDT, ω₁);

julia> result = symmetric_power(rep, 2);

julia> length(result) >= 1
true
```
"""
function symmetric_power(rep::IrrepLevi{MDT,K}, k::Integer) where {MDT,K}
  k = Int(k)
  LT = levi_type(MDT)

  k < 0 && return Vector{IrrepLevi{MDT,K}}()
  k == 0 && return [IrrepLevi(MDT, zero(rep.central),
                              WeightLatticeElem(LT === nothing ? TypeA{1} : LT,
                                               zeros(Int, rank(LT === nothing ? TypeA{1} : LT))))]
  k == 1 && return [rep]

  new_central = k * rep.central
  ss = semisimple_part(rep)

  if LT === nothing || iszero(ss)
    return [IrrepLevi(MDT, new_central, ss)]
  end

  χ = Lie.symmetric_power(ss, k)

  result = Vector{IrrepLevi{MDT,K}}()
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(MDT, new_central, hw))
    end
  end
  return result
end

# ─── Equality ────────────────────────────────────────────────────────────────

function Base.:(==)(a::IrrepLevi{MDT}, b::IrrepLevi{MDT}) where {MDT}
  p_dominant_weight(a) == p_dominant_weight(b)
end

Base.hash(a::IrrepLevi, h::UInt) = hash(a.λ, h)
