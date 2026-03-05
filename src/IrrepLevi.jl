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
    IrrepLevi{MDT}

An irreducible representation of the Levi subgroup associated to the
marked Dynkin type `MDT`, encoded by its P-dominant weight `λ` in the
weight lattice of the ambient group ``G``.

Internally the weight also stores the decomposition into:
- `central::Vector{Int}`: scaled coordinates of the central character,
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
struct IrrepLevi{MDT<:MarkedDynkinType}
  λ::WeightLatticeElem          # the P-dominant weight (ambient weight lattice)
  central::Vector{Int}           # scaled central character (×central_scaling_factor)
  semisimple::WeightLatticeElem  # highest weight of semisimple part (Levi lattice)
end

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
  C = Lie._cartan_matrix_data(DT)
  Crat = Rational{Int}.(C)
  Cinv = inv(Crat)
  sf = 1
  for j in Marked
    for k in 1:R
      sf = lcm(sf, denominator(Cinv[j, k]))
    end
  end
  return :($sf)
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
  return Rational{Int}[c // sf for c in rep.central]
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
  M = decomposition_matrix(MDT)
  Marked = marked_nodes(MDT)
  unmarked = unmarked_nodes(MDT)
  LT = levi_type(MDT)
  sf = central_scaling_factor(MDT)

  # Apply change of basis
  R = rank(_ambient_type(MDT))
  λ_vec = SVector{R,Rational{Int}}(Tuple(coefficients(λ)))
  new_coords = M * λ_vec

  # Extract central part (at marked node positions), scaled to Int
  central = Int[Int(new_coords[m] * sf) for m in Marked]

  # Extract semisimple part (at unmarked node positions)
  if LT === nothing
    # Full flag: trivial semisimple part
    semisimple = WeightLatticeElem(TypeA{1}, [0])
  else
    LR = rank(LT)
    # Natural coords: ss_nat[i] = new_coords at the i-th unmarked node
    ss_nat = [Int(new_coords[u]) for u in unmarked]
    # Apply levi_permutation: canonical coord j = ss_nat[perm[j]]
    perm = levi_permutation(MDT)
    ss_coords = [ss_nat[perm[j]] for j in 1:LR]
    semisimple = WeightLatticeElem(LT, ss_coords)
  end

  return IrrepLevi{MDT}(λ, central, semisimple)
end

"""
    IrrepLevi(::Type{MDT}, central::Vector{Int}, semisimple) -> IrrepLevi{MDT}

Construct an `IrrepLevi` from its **scaled** central character `central`
(a `Vector{Int}`, pre-multiplied by [`central_scaling_factor`](@ref)),
and the highest weight `semisimple` of its semisimple part.
The ambient P-dominant weight ``\\lambda`` is recovered automatically.
"""
function IrrepLevi(::Type{MDT}, central::Vector{Int}, semisimple::WeightLatticeElem) where {
  MDT<:MarkedDynkinType
}
  # Recover the ambient weight λ from (central, semisimple)
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  R = rank(DT)
  sf = central_scaling_factor(MDT)
  Minv = decomposition_matrix_inv(MDT)
  unmarked = unmarked_nodes(MDT)

  coords = zeros(Rational{Int}, R)
  for (idx, m) in enumerate(Marked)
    coords[m] = central[idx] // sf
  end

  # Only process semisimple part if there are unmarked nodes
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
        coords[u] = Rational{Int}(ss_vec[inv_perm[i]])
      end
    end
  end

  ambient_coords = Minv * SVector{R,Rational{Int}}(Tuple(coords))
  ambient_int = [Int(round(c)) for c in ambient_coords]
  λ = WeightLatticeElem(DT, ambient_int)

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
function tensor_product(a::IrrepLevi{MDT}, b::IrrepLevi{MDT}) where {MDT}
  LT = levi_type(MDT)

  # Central parts add (in scaled representation)
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

  result = IrrepLevi{MDT}[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(MDT, copy(new_central), hw))
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
function exterior_power(rep::IrrepLevi{MDT}, k::Integer) where {MDT}
  k = Int(k)
  LT = levi_type(MDT)

  k < 0 && return IrrepLevi{MDT}[]
  k == 0 && return [IrrepLevi(MDT, zeros(Int, length(rep.central)),
                              WeightLatticeElem(LT === nothing ? TypeA{1} : LT,
                                               zeros(Int, rank(LT === nothing ? TypeA{1} : LT))))]
  k == 1 && return [rep]

  new_central = k * rep.central
  ss = semisimple_part(rep)

  if LT === nothing || iszero(ss)
    return IrrepLevi{MDT}[]  # ⋀^k of a 1-dim rep for k > 1 is 0
  end

  dim_ss = Int(degree(ss))
  k > dim_ss && return IrrepLevi{MDT}[]

  χ = Lie.exterior_power(ss, k)

  result = IrrepLevi{MDT}[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(MDT, copy(new_central), hw))
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
function symmetric_power(rep::IrrepLevi{MDT}, k::Integer) where {MDT}
  k = Int(k)
  LT = levi_type(MDT)

  k < 0 && return IrrepLevi{MDT}[]
  k == 0 && return [IrrepLevi(MDT, zeros(Int, length(rep.central)),
                              WeightLatticeElem(LT === nothing ? TypeA{1} : LT,
                                               zeros(Int, rank(LT === nothing ? TypeA{1} : LT))))]
  k == 1 && return [rep]

  new_central = k * rep.central
  ss = semisimple_part(rep)

  if LT === nothing || iszero(ss)
    return [IrrepLevi(MDT, new_central, ss)]
  end

  χ = Lie.symmetric_power(ss, k)

  result = IrrepLevi{MDT}[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(MDT, copy(new_central), hw))
    end
  end
  return result
end

# ─── Equality ────────────────────────────────────────────────────────────────

function Base.:(==)(a::IrrepLevi{MDT}, b::IrrepLevi{MDT}) where {MDT}
  p_dominant_weight(a) == p_dominant_weight(b)
end

Base.hash(a::IrrepLevi, h::UInt) = hash(a.λ, h)
