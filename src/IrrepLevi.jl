export IrrepLevi
export central_part, semisimple_part
export to_ambient_weight, fiber_dimension, p_dominant_weight

"""
    IrrepLevi(mdt::MarkedDynkinType, λ::WeightLatticeElem)

An irreducible representation of the Levi factor attached to `mdt`, stored by
its ambient ``\\mathrm{P}``-dominant weight together with its central and semisimple Levi

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(3, 6);

julia> rep = IrrepLevi(marked_dynkin_type(X), [0, 1, -1, 0, 0]);

julia> rep == only(components(universal_subbundle(X)))
true

julia> semisimple_part(rep)
ω2

julia> fiber_dimension(rep)
3
```
coordinates.
"""
struct IrrepLevi
  mdt::MarkedDynkinType
  λ::WeightLatticeElem
  central::Vector{Int}
  semisimple::WeightLatticeElem
end

marked_dynkin_type(rep::IrrepLevi) = rep.mdt

"""
    p_dominant_weight(rep::IrrepLevi) -> WeightLatticeElem

Return the ambient ``\\mathrm{P}``-dominant weight ``\\lambda`` of the irreducible Levi
representation `rep`. This is the weight of ``\\mathrm{G}`` in fundamental weight
coordinates from which the central and semisimple parts were derived via
the decomposition matrix.
"""
p_dominant_weight(rep::IrrepLevi) = rep.λ

"""
    semisimple_part(rep::IrrepLevi) -> WeightLatticeElem

Return the highest weight of the semisimple Levi factor component of `rep`.
"""
semisimple_part(rep::IrrepLevi) = rep.semisimple

@inline function _apply_central_ext(mdt::MarkedDynkinType, λ_ivec::AbstractVector{Int})
  M = decomposition_matrix(mdt)
  sf = Int(central_scaling_factor(mdt))
  Int[
    sum(round(Int, M[m, k] * sf) * λ_ivec[k] for k in eachindex(λ_ivec)) for
    m in marked_nodes(mdt)
  ]
end

@inline function _amb_scalars(mdt::MarkedDynkinType)
  sf = Int(central_scaling_factor(mdt))
  Minv = decomposition_matrix_inv(mdt)
  sf_total = reduce(lcm, denominator.(Minv); init=sf)
  (sf_total, sf_total ÷ sf)
end

@inline function _apply_Minv_int(mdt::MarkedDynkinType, x::AbstractVector{Int})
  Minv = decomposition_matrix_inv(mdt)
  sf_total, _ = _amb_scalars(mdt)
  Int[
    sum(round(Int, Minv[i, j] * sf_total) * x[j] for j in eachindex(x)) for
    i in eachindex(x)
  ]
end

"""
    central_part(rep::IrrepLevi) -> Vector{Rational{Int}}

Return the central character of `rep` as a vector of rational numbers, one per
marked node.
"""
function central_part(rep::IrrepLevi)
  sf = central_scaling_factor(marked_dynkin_type(rep))
  Rational{Int}[c // sf for c in rep.central]
end

function _trivial_semisimple_weight(mdt::MarkedDynkinType)
  WeightLatticeElem(is_borel(mdt) ? TypeA{1} : levi_type(mdt))
end

function IrrepLevi(mdt::MarkedDynkinType, λ::WeightLatticeElem)
  unmarked = unmarked_nodes(mdt)
  λ_ivec = collect(Int, coefficients(λ))
  central = _apply_central_ext(mdt, λ_ivec)

  semisimple = if is_borel(mdt)
    _trivial_semisimple_weight(mdt)
  else
    LT = levi_type(mdt)
    perm = levi_permutation(mdt)
    WeightLatticeElem(LT, Int[λ_ivec[perm[j]] for j in 1:rank(LT)])
  end

  IrrepLevi(mdt, λ, central, semisimple)
end

function IrrepLevi(
  mdt::MarkedDynkinType, central::AbstractVector{Int}, semisimple::WeightLatticeElem
)
  length(central) == central_rank(mdt) || throw(
    ArgumentError(
      "Expected $(central_rank(mdt)) central coordinates, got $(length(central))."
    ),
  )

  DT = dynkin_type(mdt)
  marked = marked_nodes(mdt)
  unmarked = unmarked_nodes(mdt)
  R = rank(DT)

  sf_total, ratio = _amb_scalars(mdt)
  coords_full = zeros(Int, R)
  coords_full[collect(marked)] .= central .* ratio

  ss_vec = coefficients(semisimple)
  if !isempty(unmarked) && !isempty(ss_vec)
    coords_full[collect(levi_permutation(mdt))] .= ss_vec .* sf_total
  end

  λ = WeightLatticeElem(DT, div.(_apply_Minv_int(mdt, coords_full), sf_total^2))

  IrrepLevi(mdt, λ, Vector{Int}(central), semisimple)
end

function IrrepLevi(
  mdt::MarkedDynkinType, central::Vector{Rational{Int}}, semisimple::WeightLatticeElem
)
  sf = central_scaling_factor(mdt)
  central_scaled = Int.(central .* sf)
  IrrepLevi(mdt, central_scaled, semisimple)
end

"""
    to_ambient_weight(rep::IrrepLevi) -> WeightLatticeElem

Return the ambient ``\\mathrm{P}``-dominant weight of `rep` in the fundamental weight basis
of the ambient Lie algebra.
"""
to_ambient_weight(rep::IrrepLevi) = p_dominant_weight(rep)

function to_ambient_weight(mdt::MarkedDynkinType, rep::IrrepLevi)
  marked_dynkin_type(rep) == mdt ||
    throw(ArgumentError("Representation belongs to a different marked Dynkin type."))
  p_dominant_weight(rep)
end

"""
    fiber_dimension(rep::IrrepLevi) -> Int

Return the dimension of the fiber of the equivariant bundle defined by `rep`,
i.e., the dimension of the irreducible representation of the semisimple Levi
factor given by the Weyl dimension formula.
"""
function fiber_dimension(rep::IrrepLevi)::Int
  is_borel(marked_dynkin_type(rep)) && return 1
  ss = semisimple_part(rep)
  iszero(ss) && return 1
  is_dominant(ss) || return 0
  Int(degree(ss))
end

function Base.show(io::IO, rep::IrrepLevi)
  print(io, "(", sprint(show, p_dominant_weight(rep)), ")")
end

# ─── Cached tensor product ──────────────────────────────────────────────────

const _TENSOR_PRODUCT_CACHE = let b = _default_cache_budget()
  LRU{Tuple{IrrepLevi,IrrepLevi},Vector{Pair{IrrepLevi,Int}}}(;
    maxsize=_cache_maxsize(b, _DEFAULT_TENSOR_FRAC),
    by=Base.summarysize,
  )
end

"""
Canonical order for an unordered cache key: tensor products of Levi
representations are symmetric (``V ⊗ W ≅ W ⊗ V``, and every consumer treats
the decomposition as a multiset), so both orders share one cache entry.
"""
function _unordered_pair(a::IrrepLevi, b::IrrepLevi)
  hash(a) <= hash(b) ? (a, b) : (b, a)
end

function _tensor_product_terms(a::IrrepLevi, b::IrrepLevi)
  get!(_TENSOR_PRODUCT_CACHE, _unordered_pair(a, b)) do
    _tensor_product_terms_uncached(a, b)
  end
end

# The following wrappers deliberately @nospecialize on the rank-parametric
# WeightLatticeElem{DT,R}: they keep the callers from being recompiled once
# per ambient rank (same latency pattern as in Cohomology.jl).
_dual_semisimple_generic(@nospecialize(ss::WeightLatticeElem)) = dual(ss)
_tensor_product_character_generic(
  @nospecialize(a::WeightLatticeElem), @nospecialize(b::WeightLatticeElem)
) = tensor_product(a, b)
_exterior_power_character_generic(@nospecialize(ss::WeightLatticeElem), k::Int) =
  exterior_power(
    ss, k
  )
_symmetric_power_character_generic(@nospecialize(ss::WeightLatticeElem), k::Int) =
  symmetric_power(
    ss, k
  )

"""
    tensor_product(a::IrrepLevi, b::IrrepLevi) -> Vector{IrrepLevi}

Return the tensor product of two irreducible Levi representations as a list
of irreducible summands (with multiplicity).

Results are cached in a thread-safe LRU cache to avoid redundant Lie
decomposition calls.
"""
function tensor_product(a::IrrepLevi, b::IrrepLevi)
  terms = _tensor_product_terms(a, b)

  IrrepLevi[rep for (rep, mult) in terms for _ in 1:mult]
end

function _tensor_product_terms_uncached(a::IrrepLevi, b::IrrepLevi)
  marked_dynkin_type(a) == marked_dynkin_type(b) || throw(
    ArgumentError(
      "tensor_product requires Levi representations with the same marked Dynkin type."
    ),
  )

  mdt = marked_dynkin_type(a)
  new_central = a.central + b.central

  ss_a = semisimple_part(a)
  ss_b = semisimple_part(b)
  if is_borel(mdt) || (iszero(ss_a) && iszero(ss_b))
    return Pair{IrrepLevi,Int}[IrrepLevi(mdt, new_central, ss_a) => 1]
  end
  if iszero(ss_a)
    return Pair{IrrepLevi,Int}[IrrepLevi(mdt, new_central, ss_b) => 1]
  end
  if iszero(ss_b)
    return Pair{IrrepLevi,Int}[IrrepLevi(mdt, new_central, ss_a) => 1]
  end

  χ = _tensor_product_character_generic(ss_a, ss_b)
  Pair{IrrepLevi,Int}[IrrepLevi(mdt, new_central, hw) => Int(mult) for (hw, mult) in χ]
end

"""
    dual(rep::IrrepLevi) -> IrrepLevi

Return the dual (contragredient) representation of `rep`.
"""
function dual(rep::IrrepLevi)
  mdt = marked_dynkin_type(rep)
  new_central = -rep.central
  ss = semisimple_part(rep)
  if is_borel(mdt) || iszero(ss)
    return IrrepLevi(mdt, new_central, ss)
  end
  IrrepLevi(mdt, new_central, _dual_semisimple_generic(ss))
end

"""
    exterior_power(rep::IrrepLevi, k::Integer) -> Vector{IrrepLevi}

Return the `k`-th exterior power of the irreducible Levi representation `rep`
as a list of irreducible summands.
"""
function exterior_power(rep::IrrepLevi, k::Integer)
  k = Int(k)
  mdt = marked_dynkin_type(rep)

  k < 0 && return IrrepLevi[]
  k == 0 && return [
    IrrepLevi(mdt, zeros(Int, length(rep.central)), _trivial_semisimple_weight(mdt))
  ]
  k == 1 && return [rep]

  new_central = k .* rep.central
  ss = semisimple_part(rep)

  if is_borel(mdt) || iszero(ss)
    return IrrepLevi[]
  end

  dim_ss = Int(degree(ss))
  k > dim_ss && return IrrepLevi[]

  χ = _exterior_power_character_generic(ss, k)
  IrrepLevi[IrrepLevi(mdt, new_central, hw) for (hw, mult) in χ for _ in 1:mult]
end

"""
    symmetric_power(rep::IrrepLevi, k::Integer) -> Vector{IrrepLevi}

Return the `k`-th symmetric power of the irreducible Levi representation `rep`
as a list of irreducible summands.
"""
function symmetric_power(rep::IrrepLevi, k::Integer)
  k = Int(k)
  mdt = marked_dynkin_type(rep)

  k < 0 && return IrrepLevi[]
  k == 0 && return [
    IrrepLevi(mdt, zeros(Int, length(rep.central)), _trivial_semisimple_weight(mdt))
  ]
  k == 1 && return [rep]

  new_central = k .* rep.central
  ss = semisimple_part(rep)
  if is_borel(mdt) || iszero(ss)
    return [IrrepLevi(mdt, new_central, ss)]
  end

  χ = _symmetric_power_character_generic(ss, k)
  IrrepLevi[IrrepLevi(mdt, new_central, hw) for (hw, mult) in χ for _ in 1:mult]
end

"""
    IrrepLevi(mdt::MarkedDynkinType, coeffs::AbstractVector{<:Integer})

Convenience constructor: build an irreducible Levi representation from a
vector of fundamental-weight coefficients ``(\\lambda_1, \\ldots, \\lambda_r)``
in the ambient weight lattice.
"""
function IrrepLevi(mdt::MarkedDynkinType, coeffs::AbstractVector{<:Integer})
  IrrepLevi(mdt, WeightLatticeElem(dynkin_type(mdt), Vector{Int}(coeffs)))
end

Base.:(==)(a::IrrepLevi, b::IrrepLevi) =
  marked_dynkin_type(a) == marked_dynkin_type(b) &&
  p_dominant_weight(a) == p_dominant_weight(b)

Base.hash(a::IrrepLevi, h::UInt) = hash(a.λ, hash(a.mdt, h))
