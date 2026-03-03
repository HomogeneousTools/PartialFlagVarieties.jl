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
- `central::Vector{Rational{Int}}`: coordinates of the central character,
  indexed by the marked (nonparabolic) nodes
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
  central::Vector{Rational{Int}} # central character coordinates (at marked nodes)
  semisimple::WeightLatticeElem  # highest weight of semisimple part (Levi lattice)
end

# ─── Core accessors ──────────────────────────────────────────────────────────

"""
    p_dominant_weight(rep::IrrepLevi) -> WeightLatticeElem

Return the P-dominant weight ``\\lambda`` that defines this Levi representation.
"""
p_dominant_weight(rep::IrrepLevi) = rep.λ

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
    IrrepLevi(::Type{MDT}, central, semisimple) -> IrrepLevi{MDT}

Construct an `IrrepLevi` from its central character `central` (a
`Vector{Rational{Int}}` indexed by marked nodes) and the highest weight
`semisimple` of its semisimple part (a `WeightLatticeElem` in the Levi
weight lattice).  The ambient P-dominant weight ``\\lambda`` is recovered
automatically.
"""
function IrrepLevi(::Type{MDT}, central::Vector{Rational{Int}}, semisimple::WeightLatticeElem) where {
  MDT<:MarkedDynkinType
}
  # Recover the ambient weight λ from (central, semisimple)
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  R = rank(DT)
  Minv = decomposition_matrix_inv(MDT)
  unmarked = unmarked_nodes(MDT)

  coords = zeros(Rational{Int}, R)
  for (idx, m) in enumerate(Marked)
    coords[m] = central[idx]
  end

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

  ambient_coords = Minv * SVector{R,Rational{Int}}(Tuple(coords))
  ambient_int = [Int(round(c)) for c in ambient_coords]
  λ = WeightLatticeElem(DT, ambient_int)

  return IrrepLevi{MDT}(λ, central, semisimple)
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
  Minv = decomposition_matrix_inv(MDT)
  unmarked = unmarked_nodes(MDT)

  # Reconstruct the full coordinate vector
  coords = zeros(Rational{Int}, R)

  # Central part at marked positions
  for (idx, m) in enumerate(Marked)
    coords[m] = rep.central[idx]
  end

  # Semisimple part at unmarked positions
  # Apply inverse of levi_permutation: ss_nat[perm[j]] = ss_canonical[j]
  ss_vec = coefficients(rep.semisimple)
  LR = length(ss_vec)
  if LR > 0
    perm = levi_permutation(MDT)
    # Invert the permutation: nat[perm[j]] = canon[j]  so  nat[i] = canon[perm⁻¹[i]]
    inv_perm = Vector{Int}(undef, LR)
    for j in 1:LR
      inv_perm[perm[j]] = j
    end
    for (i, u) in enumerate(unmarked)
      coords[u] = Rational{Int}(ss_vec[inv_perm[i]])
    end
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
2

julia> ω₂ = fundamental_weight(TypeA{3}, 2);

julia> fiber_dimension(IrrepLevi(MDT, ω₂))
2
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
  # Print the P-dominant weight directly
  print(io, "(", sprint(show, rep.λ), ")")
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
    return [IrrepLevi(MDT, new_central, a.semisimple)]
  end

  if iszero(a.semisimple)
    return [IrrepLevi(MDT, new_central, b.semisimple)]
  end
  if iszero(b.semisimple)
    return [IrrepLevi(MDT, new_central, a.semisimple)]
  end

  # Tensor product via Lie.jl's Weyl character ring
  χ = Lie.tensor_product(a.semisimple, b.semisimple)

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

julia> d.central == -rep.central
true
```
"""
function dual(rep::IrrepLevi{MDT}) where {MDT}
  LT = levi_type(MDT)
  new_central = -rep.central

  if LT === nothing || iszero(rep.semisimple)
    return IrrepLevi(MDT, new_central, rep.semisimple)
  end

  new_ss = Lie.dual(rep.semisimple)
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
function exterior_power(rep::IrrepLevi{MDT}, k::Int) where {MDT}
  LT = levi_type(MDT)

  k < 0 && return IrrepLevi{MDT}[]
  k == 0 && return [IrrepLevi(MDT, zeros(Rational{Int}, length(rep.central)),
                              WeightLatticeElem(LT === nothing ? TypeA{1} : LT,
                                               zeros(Int, rank(LT === nothing ? TypeA{1} : LT))))]
  k == 1 && return [rep]

  new_central = k * rep.central

  if LT === nothing || iszero(rep.semisimple)
    return IrrepLevi{MDT}[]  # ⋀^k of a 1-dim rep for k > 1 is 0
  end

  dim_ss = Int(degree(rep.semisimple))
  k > dim_ss && return IrrepLevi{MDT}[]

  χ = Lie.exterior_power(rep.semisimple, k)

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
function symmetric_power(rep::IrrepLevi{MDT}, k::Int) where {MDT}
  LT = levi_type(MDT)

  k < 0 && return IrrepLevi{MDT}[]
  k == 0 && return [IrrepLevi(MDT, zeros(Rational{Int}, length(rep.central)),
                              WeightLatticeElem(LT === nothing ? TypeA{1} : LT,
                                               zeros(Int, rank(LT === nothing ? TypeA{1} : LT))))]
  k == 1 && return [rep]

  new_central = k * rep.central

  if LT === nothing || iszero(rep.semisimple)
    return [IrrepLevi(MDT, new_central, rep.semisimple)]
  end

  χ = Lie.symmetric_power(rep.semisimple, k)

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
  return a.λ == b.λ
end

Base.hash(a::IrrepLevi, h::UInt) = hash(a.λ, h)
