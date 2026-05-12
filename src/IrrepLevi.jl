export IrrepLevi
export central_part, semisimple_part
export to_ambient_weight, fiber_dimension, p_dominant_weight

"""
    IrrepLevi(mdt::MarkedDynkinType, λ::WeightLatticeElem)

An irreducible representation of the Levi factor attached to `mdt`, stored by
its ambient ``P``-dominant weight together with its central and semisimple Levi

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

Return the ambient ``P``-dominant weight ``\\lambda`` of the irreducible Levi
representation `rep`. This is the weight of ``G`` in fundamental weight
coordinates from which the central and semisimple parts were derived via
the decomposition matrix.
"""
p_dominant_weight(rep::IrrepLevi) = rep.λ

"""
    semisimple_part(rep::IrrepLevi) -> WeightLatticeElem

Return the highest weight of the semisimple Levi factor component of `rep`.
"""
semisimple_part(rep::IrrepLevi) = rep.semisimple

_central_length(mdt::MarkedDynkinType) = central_rank(mdt)

@inline function _apply_central_ext(mdt::MarkedDynkinType, λ_ivec::AbstractVector{Int})
  marked = marked_nodes(mdt)
  M = decomposition_matrix(mdt)
  sf = Int(central_scaling_factor(mdt))
  central = Vector{Int}(undef, length(marked))
  for (idx, m) in enumerate(marked)
    total = 0
    for k in eachindex(λ_ivec)
      total += round(Int, M[m, k] * sf) * λ_ivec[k]
    end
    central[idx] = total
  end
  central
end

@inline function _amb_scalars(mdt::MarkedDynkinType)
  sf = Int(central_scaling_factor(mdt))
  Minv = decomposition_matrix_inv(mdt)
  sf_total = sf
  for j in axes(Minv, 1), k in axes(Minv, 2)
    sf_total = lcm(sf_total, denominator(Minv[j, k]))
  end
  (sf_total, sf_total ÷ sf)
end

@inline function _apply_Minv_int(mdt::MarkedDynkinType, x::AbstractVector{Int})
  Minv = decomposition_matrix_inv(mdt)
  sf_total, _ = _amb_scalars(mdt)
  result = Vector{Int}(undef, length(x))
  for i in eachindex(result)
    total = 0
    for j in eachindex(x)
      total += round(Int, Minv[i, j] * sf_total) * x[j]
    end
    result[i] = total
  end
  result
end

"""
    central_part(rep::IrrepLevi) -> Vector{Rational{Int}}

Return the central character of `rep` as a vector of rational numbers, one per
marked node.
"""
function central_part(rep::IrrepLevi)
  sf = central_scaling_factor(marked_dynkin_type(rep))
  Vector{Rational{Int}}(Rational{Int}[c // sf for c in rep.central])
end

function _trivial_semisimple_weight(mdt::MarkedDynkinType)
  trivial_type = is_borel(mdt) ? TypeA{1} : levi_type(mdt)
  WeightLatticeElem(trivial_type, zeros(Int, rank(trivial_type)))
end

function IrrepLevi(mdt::MarkedDynkinType, λ::WeightLatticeElem)
  unmarked = unmarked_nodes(mdt)
  λ_coeffs = coefficients(λ)
  λ_ivec = Vector{Int}(undef, length(λ_coeffs))
  for i in eachindex(λ_coeffs)
    λ_ivec[i] = λ_coeffs[i]
  end
  central = _apply_central_ext(mdt, λ_ivec)

  if is_borel(mdt)
    semisimple = _trivial_semisimple_weight(mdt)
  else
    LT = levi_type(mdt)
    LR = rank(LT)
    perm = levi_permutation(mdt)
    ss_coords = Vector{Int}(undef, LR)
    for j in 1:LR
      ss_coords[j] = λ_ivec[unmarked[perm[j]]]
    end
    semisimple = WeightLatticeElem(LT, ss_coords)
  end

  IrrepLevi(mdt, λ, central, semisimple)
end

function IrrepLevi(
  mdt::MarkedDynkinType, central::AbstractVector{Int}, semisimple::WeightLatticeElem
)
  length(central) == _central_length(mdt) || throw(
    ArgumentError(
      "Expected $(_central_length(mdt)) central coordinates, got $(length(central))."
    ),
  )

  DT = dynkin_type(mdt)
  marked = marked_nodes(mdt)
  unmarked = unmarked_nodes(mdt)
  R = rank(DT)

  sf_total, ratio = _amb_scalars(mdt)
  coords_full = Vector{Int}(undef, R)

  for (idx, m) in enumerate(marked)
    coords_full[m] = Int(central[idx]) * ratio
  end

  if !isempty(unmarked)
    ss_vec = coefficients(semisimple)
    LR = length(ss_vec)
    if LR > 0
      perm = levi_permutation(mdt)
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

  ambient_scaled = _apply_Minv_int(mdt, coords_full)
  sf_sq = sf_total * sf_total
  λ_coords = Vector{Int}(undef, length(ambient_scaled))
  for i in eachindex(ambient_scaled)
    λ_coords[i] = div(ambient_scaled[i], sf_sq)
  end
  λ = WeightLatticeElem(DT, λ_coords)

  IrrepLevi(mdt, λ, Vector{Int}(central), semisimple)
end

function IrrepLevi(
  mdt::MarkedDynkinType, central::Vector{Rational{Int}}, semisimple::WeightLatticeElem
)
  sf = central_scaling_factor(mdt)
  central_scaled = Int[Int(c * sf) for c in central]
  IrrepLevi(mdt, central_scaled, semisimple)
end

"""
    to_ambient_weight(rep::IrrepLevi) -> WeightLatticeElem

Return the ambient `P`-dominant weight of `rep` in the fundamental weight basis
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

function _tensor_product_terms(a::IrrepLevi, b::IrrepLevi)
  get!(_TENSOR_PRODUCT_CACHE, (a, b)) do
    _tensor_product_terms_uncached(a, b)
  end
end

_dual_semisimple_generic(@nospecialize(ss::WeightLatticeElem)) = dual(ss)
_tensor_product_character_generic(
@nospecialize(a::WeightLatticeElem), @nospecialize(b::WeightLatticeElem)
) = tensor_product(a, b)
_exterior_power_character_generic(@nospecialize(ss::WeightLatticeElem), k::Int) = exterior_power(
  ss, k
)
_symmetric_power_character_generic(@nospecialize(ss::WeightLatticeElem), k::Int) = symmetric_power(
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

  result = IrrepLevi[]
  for (rep, mult) in terms
    for _ in 1:mult
      push!(result, rep)
    end
  end
  result
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
  result = Pair{IrrepLevi,Int}[]
  for (hw, mult) in χ
    push!(result, IrrepLevi(mdt, new_central, hw) => Int(mult))
  end
  result
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
  result = IrrepLevi[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(mdt, new_central, hw))
    end
  end
  result
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
  result = IrrepLevi[]
  for (hw, mult) in χ
    for _ in 1:mult
      push!(result, IrrepLevi(mdt, new_central, hw))
    end
  end
  result
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
