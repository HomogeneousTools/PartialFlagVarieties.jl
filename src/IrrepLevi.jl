# ═══════════════════════════════════════════════════════════════════════════════
#  IrrepLevi — irreducible representation of the Levi subgroup
#
#  An irreducible representation of the Levi subgroup L of a parabolic P
#  decomposes as a character of the center Z(L) (the "central part")
#  tensored with an irreducible representation of the semisimple part [L,L]
#  (the "semisimple part").
#
#  The decomposition uses the special_matrix change-of-basis from
#  MarkedDynkinType.jl.
# ═══════════════════════════════════════════════════════════════════════════════

export IrrepLevi
export central_part, semisimple_part
export to_ambient_weight, fiber_dimension

# Names from Lie and StaticArrays are available via the parent module's
# `using Lie` and `using StaticArrays`.

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    IrrepLevi{MDT}

An irreducible representation of the Levi subgroup associated to the
marked Dynkin type `MDT`.

# Fields
- `central::Vector{Rational{Int}}`: coordinates of the central character,
  indexed by the marked (nonparabolic) nodes
- `semisimple::WeightLatticeElem`: highest weight of the semisimple part,
  as a weight in the Levi's weight lattice

The ambient weight ``\\lambda`` decomposes under the Levi as:
``\\lambda \\mapsto (\\text{central part}, \\text{semisimple part})``
via the [`special_matrix`](@ref) change of basis.

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
  central::Vector{Rational{Int}}
  semisimple::WeightLatticeElem
end

# ─── Core accessors ──────────────────────────────────────────────────────────

"""
    central_part(rep::IrrepLevi) -> Vector{Rational{Int}}

Return the central character part of the Levi representation.
"""
central_part(rep::IrrepLevi) = rep.central

"""
    semisimple_part(rep::IrrepLevi) -> WeightLatticeElem

Return the highest weight of the semisimple part.
"""
semisimple_part(rep::IrrepLevi) = rep.semisimple

# ─── Construction from ambient weight ───────────────────────────────────────

"""
    IrrepLevi(::Type{MDT}, λ::WeightLatticeElem) -> IrrepLevi{MDT}

Construct an irreducible Levi representation from a weight ``\\lambda``
of the ambient group ``G``, by applying the [`special_matrix`](@ref)
change of basis to decompose into central + semisimple parts.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (2,)};

julia> ω₁ = fundamental_weight(TypeA{4}, 1);

julia> rep = IrrepLevi(MDT, ω₁);

julia> fiber_dimension(rep)
1
```
"""
function IrrepLevi(::Type{MDT}, λ::WeightLatticeElem) where {
  MDT<:MarkedDynkinType
}
  M = special_matrix(MDT)
  Marked = marked_nodes(MDT)
  um = unmarked_nodes(MDT)
  LT = levi_type(MDT)

  # Apply change of basis
  R = rank(_ambient_type(MDT))
  λ_vec = SVector{R,Rational{Int}}(Tuple(coefficients(λ)))
  new_coords = M * λ_vec

  # Extract central part (at marked node positions)
  central = Rational{Int}[new_coords[m] for m in Marked]

  # Extract semisimple part (at unmarked node positions)
  if LT === nothing
    # Full flag: trivial semisimple part
    semisimple = WeightLatticeElem(TypeA{1}, [0])
  else
    LR = rank(LT)
    ss_coords = [Int(new_coords[u]) for u in um]
    semisimple = WeightLatticeElem(LT, ss_coords)
  end

  return IrrepLevi{MDT}(central, semisimple)
end

# ─── Back-conversion to ambient weight ───────────────────────────────────────

"""
    to_ambient_weight(::Type{MDT}, rep::IrrepLevi{MDT}) -> WeightLatticeElem

Convert an `IrrepLevi` back to a weight in the ambient group's weight lattice.

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
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  R = rank(DT)
  Minv = special_matrix_inv(MDT)
  um = unmarked_nodes(MDT)

  # Reconstruct the full coordinate vector
  coords = zeros(Rational{Int}, R)

  # Central part at marked positions
  for (idx, m) in enumerate(Marked)
    coords[m] = rep.central[idx]
  end

  # Semisimple part at unmarked positions
  ss_vec = coefficients(rep.semisimple)
  for (idx, u) in enumerate(um)
    coords[u] = Rational{Int}(ss_vec[idx])
  end

  # Apply inverse change of basis
  ambient_coords = Minv * SVector{R,Rational{Int}}(Tuple(coords))

  # Convert to Int (should be integral)
  ambient_int = [Int(round(c)) for c in ambient_coords]

  return WeightLatticeElem(DT, ambient_int)
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
1

julia> ω₂ = fundamental_weight(TypeA{3}, 2);

julia> fiber_dimension(IrrepLevi(MDT, ω₂))
1
```
"""
function fiber_dimension(rep::IrrepLevi{MDT}) where {MDT}
  LT = levi_type(MDT)
  LT === nothing && return BigInt(1)
  iszero(rep.semisimple) && return BigInt(1)
  is_dominant(rep.semisimple) || return BigInt(0)
  return degree(rep.semisimple)
end

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, rep::IrrepLevi{MDT}) where {MDT}
  # Format central part
  central_strs = String[]
  for c in rep.central
    if denominator(c) == 1
      push!(central_strs, string(numerator(c)))
    else
      push!(central_strs, string(c))
    end
  end
  central_str = join(central_strs, ", ")

  # Format semisimple part
  ss_str = sprint(show, rep.semisimple)

  print(io, "(", central_str, " | ", ss_str, ")")
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

  # Central parts add
  new_central = a.central + b.central

  # Semisimple parts: tensor product decomposition
  if LT === nothing || (iszero(a.semisimple) && iszero(b.semisimple))
    return [IrrepLevi{MDT}(new_central, a.semisimple)]
  end

  if iszero(a.semisimple)
    return [IrrepLevi{MDT}(new_central, b.semisimple)]
  end
  if iszero(b.semisimple)
    return [IrrepLevi{MDT}(new_central, a.semisimple)]
  end

  # Tensor product via Lie.jl's Weyl character ring
  χ = Lie.tensor_product(a.semisimple, b.semisimple)

  result = IrrepLevi{MDT}[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi{MDT}(copy(new_central), hw))
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

julia> d.central == -rep.central
true
```
"""
function dual(rep::IrrepLevi{MDT}) where {MDT}
  LT = levi_type(MDT)
  new_central = -rep.central

  if LT === nothing || iszero(rep.semisimple)
    return IrrepLevi{MDT}(new_central, rep.semisimple)
  end

  new_ss = Lie.dual(rep.semisimple)
  return IrrepLevi{MDT}(new_central, new_ss)
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
function exterior_power(rep::IrrepLevi{MDT}, k::Int) where {MDT}
  LT = levi_type(MDT)

  k < 0 && return IrrepLevi{MDT}[]
  if k == 0
    zero_central = zeros(Rational{Int}, length(rep.central))
    LT_actual = LT === nothing ? TypeA{1} : LT
    zero_ss = WeightLatticeElem(LT_actual, zeros(Int, rank(LT_actual)))
    return [IrrepLevi{MDT}(zero_central, zero_ss)]
  end
  k == 1 && return [IrrepLevi{MDT}(copy(rep.central), rep.semisimple)]

  new_central = k * rep.central

  if LT === nothing || iszero(rep.semisimple)
    if k == 1
      return [IrrepLevi{MDT}(new_central, rep.semisimple)]
    else
      return IrrepLevi{MDT}[]  # ⋀^k of a 1-dim rep for k > 1 is 0
    end
  end

  dim_ss = Int(degree(rep.semisimple))
  k > dim_ss && return IrrepLevi{MDT}[]

  χ = Lie.exterior_power(rep.semisimple, k)

  result = IrrepLevi{MDT}[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi{MDT}(copy(new_central), hw))
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
function symmetric_power(rep::IrrepLevi{MDT}, k::Int) where {MDT}
  LT = levi_type(MDT)

  k < 0 && return IrrepLevi{MDT}[]
  if k == 0
    zero_central = zeros(Rational{Int}, length(rep.central))
    LT_actual = LT === nothing ? TypeA{1} : LT
    zero_ss = WeightLatticeElem(LT_actual, zeros(Int, rank(LT_actual)))
    return [IrrepLevi{MDT}(zero_central, zero_ss)]
  end
  k == 1 && return [IrrepLevi{MDT}(copy(rep.central), rep.semisimple)]

  new_central = k * rep.central

  if LT === nothing || iszero(rep.semisimple)
    return [IrrepLevi{MDT}(new_central, rep.semisimple)]
  end

  χ = Lie.symmetric_power(rep.semisimple, k)

  result = IrrepLevi{MDT}[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi{MDT}(copy(new_central), hw))
    end
  end
  return result
end

# ─── Equality ────────────────────────────────────────────────────────────────

function Base.:(==)(a::IrrepLevi{MDT}, b::IrrepLevi{MDT}) where {MDT}
  return a.central == b.central && a.semisimple == b.semisimple
end

Base.hash(a::IrrepLevi, h::UInt) = hash((a.central, a.semisimple), h)
