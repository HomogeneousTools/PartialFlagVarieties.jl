# ═══════════════════════════════════════════════════════════════════════════════
#  CompletelyReducibleBundle — a completely reducible equivariant bundle on G/P
#
#  A completely reducible equivariant bundle on G/P is an equivariant bundle
#  whose underlying P-representation is completely reducible (semisimple).
#  It is encoded as a formal (virtual) sum of irreducible Levi representations,
#  together with a reference to the underlying partial flag variety.
# ═══════════════════════════════════════════════════════════════════════════════

export Bundle
export CompletelyReducibleBundle
export components, variety
export rank_bundle, tangent_bundle, cotangent_bundle
export structure_sheaf, O, zero_bundle, line_bundle, canonical_bundle, anticanonical_bundle
export T, E  # shorthands for tangent_bundle and the CompletelyReducibleBundle constructors
export det, determinant
export picard_degrees
export fano_index

# Names from Lie, StaticArrays, Combinatorics are available via the parent module.

# ═══════════════════════════════════════════════════════════════════════════════
#  Abstract Bundle type
# ═══════════════════════════════════════════════════════════════════════════════

"""
  Bundle

Abstract supertype for equivariant vector bundles on a partial flag variety.

Concrete subtypes:
- [`CompletelyReducibleBundle`](@ref): semisimple equivariant bundles
- [`FilteredBundle`](@ref): bundles with a filtration by equivariant subbundles
"""
abstract type Bundle end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
  CompletelyReducibleBundle

A completely reducible equivariant vector bundle on a partial flag variety.

Stored as a list of irreducible Levi representations (with multiplicity
encoded by repetition), together with the underlying [`PartialFlagVariety`](@ref).

# Fields
- `variety::PartialFlagVariety`: the partial flag variety
- `components::Vector{IrrepLevi}`: the irreducible summands

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> T = tangent_bundle(X);

julia> rank_bundle(T)
4
```
"""
struct CompletelyReducibleBundle <: Bundle
  variety::PartialFlagVariety
  components::Vector{IrrepLevi}

  function CompletelyReducibleBundle(
    variety::PartialFlagVariety,
    components::Vector{IrrepLevi},
  )
    bundle_mdt = marked_dynkin_type(variety)
    for (idx, component) in enumerate(components)
      marked_dynkin_type(component) == bundle_mdt || throw(
        ArgumentError(
          "Bundle component $idx belongs to $(marked_dynkin_type(component)), expected $bundle_mdt."
        ),
      )
    end
    new(variety, components)
  end
end

Base.:(==)(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle) =
  E.variety == F.variety && E.components == F.components
Base.hash(E::CompletelyReducibleBundle, h::UInt) = hash(E.components, hash(E.variety, h))

"""
    CompletelyReducibleBundle(X::PartialFlagVariety, λ::WeightLatticeElem) -> CompletelyReducibleBundle

Construct the equivariant bundle on the partial flag variety `X``
corresponding to the Levi representation of highest weight `λ`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(3, 6);

julia> λ = WeightLatticeElem(dynkin_type(X), [0, 1, -1, 0, 0]);

julia> E = CompletelyReducibleBundle(X, λ);

julia> components(E) == components(universal_subbundle(X))
true
```
"""
function CompletelyReducibleBundle(X::PartialFlagVariety, λ::WeightLatticeElem)
  CompletelyReducibleBundle(X, [λ])
end

"""
    CompletelyReducibleBundle(X::PartialFlagVariety, weights::AbstractVector{<:WeightLatticeElem})

Convenience constructor: build the direct sum of the equivariant bundles on
`X` corresponding to the ambient highest weights in `weights`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(3, 6);

julia> weights = [WeightLatticeElem(dynkin_type(X), [0, 1, -1, 0, 0]), WeightLatticeElem(dynkin_type(X))];

julia> E = CompletelyReducibleBundle(X, weights);

julia> components(E) == vcat(components(universal_subbundle(X)), components(structure_sheaf(X)))
true
```
"""
function CompletelyReducibleBundle(
  X::PartialFlagVariety, weights::AbstractVector{<:WeightLatticeElem}
)
  mdt = marked_dynkin_type(X)
  components = [IrrepLevi(mdt, λ) for λ in weights]
  CompletelyReducibleBundle(X, components)
end

"""
    CompletelyReducibleBundle(X::PartialFlagVariety, coeffs::AbstractVector{<:Integer})

Convenience constructor: build the equivariant bundle on `X` corresponding to
the Levi representation whose ambient highest weight has the given
fundamental-weight coefficients.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(3, 6);

julia> U = CompletelyReducibleBundle(X, [0, 1, -1, 0, 0]);

julia> components(U) == components(universal_subbundle(X))
true

julia> rank_bundle(U)
3
```
"""
function CompletelyReducibleBundle(X::PartialFlagVariety, coeffs::AbstractVector{<:Integer})
  CompletelyReducibleBundle(X, WeightLatticeElem(dynkin_type(X), Vector{Int}(coeffs)))
end

"""
    CompletelyReducibleBundle(X::PartialFlagVariety, coeffs_list::AbstractVector{<:AbstractVector{<:Integer}})

Convenience constructor: build the direct sum of the equivariant bundles on
`X` corresponding to the ambient highest weights given by the vectors of
fundamental-weight coefficients in `coeffs_list`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(3, 6);

julia> E = CompletelyReducibleBundle(X, [[0, 1, -1, 0, 0], [0, 0, 0, 0, 0]]);

julia> components(E) == vcat(components(universal_subbundle(X)), components(structure_sheaf(X)))
true

julia> rank_bundle(E)
4
```
"""
function CompletelyReducibleBundle(
  X::PartialFlagVariety, coeffs_list::AbstractVector{<:AbstractVector{<:Integer}}
)
  mdt = marked_dynkin_type(X)
  components = [IrrepLevi(mdt, coeffs) for coeffs in coeffs_list]
  CompletelyReducibleBundle(X, components)
end

"""
    E(X::PartialFlagVariety, weights) -> CompletelyReducibleBundle

Shorthand for the [`CompletelyReducibleBundle`](@ref) constructors: build the
equivariant bundle on `X` from `weights`, which may be a single
`WeightLatticeElem`, a vector of them, a vector of fundamental-weight
coefficients, or a list of such vectors.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(3, 6);

julia> components(E(X, [0, 1, -1, 0, 0])) == components(universal_subbundle(X))
true
```
"""
E(X::PartialFlagVariety, weights) = CompletelyReducibleBundle(X, weights)

# ─── Accessors ───────────────────────────────────────────────────────────────

"""
    variety(E::CompletelyReducibleBundle) -> PartialFlagVariety

Return the partial flag variety on which this bundle lives.
"""
variety(E::CompletelyReducibleBundle) = E.variety

"""
    components(E::CompletelyReducibleBundle) -> Vector{<:IrrepLevi}

Return the irreducible summands.
"""
components(E::CompletelyReducibleBundle) = E.components

"""
    n_components(E::CompletelyReducibleBundle) -> Int

Return the number of irreducible summands.
"""
n_components(E::CompletelyReducibleBundle) = length(E.components)

"""
    rank_bundle(E::CompletelyReducibleBundle) -> Int

Return the total rank (fiber dimension) of the bundle.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> rank_bundle(structure_sheaf(X))
1
```
"""
function rank_bundle(E::CompletelyReducibleBundle)
  sum(fiber_dimension(component) for component in E.components; init=0)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Standard bundles
# ═══════════════════════════════════════════════════════════════════════════════

"""
    O(X::PartialFlagVariety) -> CompletelyReducibleBundle
    structure_sheaf(X::PartialFlagVariety) -> CompletelyReducibleBundle

The trivial line bundle ``\\mathcal{O}_{G/P}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> rank_bundle(structure_sheaf(X))
1
```
"""
function O(X::PartialFlagVariety)
  CompletelyReducibleBundle(X, WeightLatticeElem(dynkin_type(X)))
end

"""
    O(X::PartialFlagVariety, d::Integer) -> CompletelyReducibleBundle

The line bundle ``\\mathcal{O}(d)`` on `X`. Shorthand for
[`line_bundle(X, d)`](@ref); requires Picard rank 1.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(O(Gr(2, 4), 3))
1
```
"""
O(X::PartialFlagVariety, d::Integer) = line_bundle(X, d)

"""
    O(X::PartialFlagVariety, degrees::Vector{<:Integer}) -> CompletelyReducibleBundle

The line bundle ``\\mathcal{O}(d_1, \\ldots, d_r)`` on `X`, one degree per marked
node. Shorthand for [`line_bundle(X, degrees)`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{3}, (1, 3));

julia> rank_bundle(O(X, [2, 1]))
1
```
"""
O(X::PartialFlagVariety, degrees::Vector{<:Integer}) = line_bundle(X, degrees)

"""
    structure_sheaf(X::PartialFlagVariety) -> CompletelyReducibleBundle

The trivial line bundle ``\\mathcal{O}_{G/P}``. Alias for `O`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> rank_bundle(structure_sheaf(X))
1
```
"""
structure_sheaf(X::PartialFlagVariety) = O(X)

"""
    zero_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The zero bundle on `X` (the empty direct sum).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> E = zero_bundle(X);

julia> rank_bundle(E)
0

julia> n_components(E)
0
```
"""
function zero_bundle(X::PartialFlagVariety)
  CompletelyReducibleBundle(X, IrrepLevi[])
end

"""
    line_bundle(X::PartialFlagVariety, i::Integer) -> CompletelyReducibleBundle

The line bundle ``\\mathcal{O}(i)`` on `X`.

When the Picard rank of `X` is 1 (generalized Grassmannian), `i` is the
degree. When the Picard rank is greater than 1, this method raises an
error — use [`line_bundle(X, degrees)`](@ref) instead.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> L = line_bundle(X, 1);

julia> rank_bundle(L)
1
```
"""
function line_bundle(X::PartialFlagVariety, i::Integer)
  picard_rank(X) == 1 || throw(
    ArgumentError(
      "line_bundle(X, i::Integer) requires Picard rank 1, but X has Picard rank $(picard_rank(X)). " *
      "Use line_bundle(X, degrees::Vector{<:Integer}) instead.",
    ),
  )

  line_bundle(X, [i])
end

"""
    line_bundle(X::PartialFlagVariety, degrees::Vector{<:Integer}) -> CompletelyReducibleBundle

The line bundle ``\\mathcal{O}(d_1, \\ldots, d_r)`` on `X`, where ``d_j`` is the
degree at the ``j``-th marked node.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{3}, (1, 3));

julia> L = line_bundle(X, [2, 1]);

julia> rank_bundle(L)
1
```
"""
function line_bundle(X::PartialFlagVariety, degrees::Vector{<:Integer})
  marked = marked_nodes(X)
  length(degrees) == length(marked) || throw(
    ArgumentError(
      "Expected $(length(marked)) degrees (one per marked node), got $(length(degrees))."
    ),
  )

  coeffs = zeros(Int, rank(X))
  coeffs[collect(marked)] .= degrees

  λ = WeightLatticeElem(dynkin_type(X), coeffs)
  CompletelyReducibleBundle(X, λ)
end

"""
    picard_degrees(E::CompletelyReducibleBundle) -> Vector{Int}

Return the Picard-group coordinates of the rank-1 bundle `E`: the
fundamental-weight degrees at each marked node, i.e. the vector `d` such
that `E == line_bundle(variety(E), d)`.

Throws an `ArgumentError` if `E` does not have rank 1.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> picard_degrees(anticanonical_bundle(X))
1-element Vector{Int64}:
 4

julia> X2 = partial_flag_variety(TypeA{3}, (1, 3));

julia> picard_degrees(line_bundle(X2, [3, 2]))
2-element Vector{Int64}:
 3
 2
```
"""
function picard_degrees(E::CompletelyReducibleBundle)
  rank_bundle(E) == 1 || throw(
    ArgumentError("picard_degrees requires a rank-1 bundle, got rank $(rank_bundle(E)).")
  )
  X = variety(E)
  v = coefficients(p_dominant_weight(only(components(E))))
  Int[v[m] for m in marked_nodes(X)]
end

"""
    is_ample_line_bundle(E::CompletelyReducibleBundle) -> Bool

Whether `E` is an ample line bundle: rank 1 with strictly positive degree
at every marked node.  Returns `false` for bundles of rank different
from 1.  The restriction of an ample line bundle to a subvariety is again
ample.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> is_ample_line_bundle(line_bundle(X, 1))
true

julia> is_ample_line_bundle(structure_sheaf(X))
false

julia> is_ample_line_bundle(universal_subbundle(X))
false
```
"""
function is_ample_line_bundle(E::CompletelyReducibleBundle)
  rank_bundle(E) == 1 && all(>(0), picard_degrees(E))
end

"""
    _product_factor_range(total_rank, factor_rank, offset) -> UnitRange{Int}

Return the block of ambient coordinates occupied by a factor embedded into a
product ambient of rank `total_rank`.
"""
function _product_factor_range(total_rank::Int, factor_rank::Int, offset::Int)
  0 <= offset <= total_rank || throw(ArgumentError("Offset $offset is out of range."))

  last_index = offset + factor_rank
  last_index <= total_rank || throw(
    ArgumentError(
      "A factor of rank $factor_rank does not fit in an ambient rank $total_rank at offset $offset."
    ),
  )
  (offset + 1):last_index
end

"""
    _embed_ambient_weight(DT, λ, offset) -> WeightLatticeElem

Embed the weight `λ` into the ambient weight lattice of `DT` by placing its
coefficients into the coordinate block determined by `offset`.
"""
function _embed_ambient_weight(
  ::Type{DT}, λ::WeightLatticeElem, offset::Int
) where {DT<:DynkinType}
  weight_coefficients = collect(Int, coefficients(λ))
  coordinate_range = _product_factor_range(rank(DT), length(weight_coefficients), offset)
  ambient_coefficients = zeros(Int, rank(DT))
  ambient_coefficients[coordinate_range] = weight_coefficients
  WeightLatticeElem(DT, ambient_coefficients)
end

"""
    _lift_irrep_to_product(X, rep, offset) -> IrrepLevi

Lift an irreducible Levi representation from a factor ambient into the product
ambient `X` by embedding its ambient highest weight at the specified offset.
"""
function _lift_irrep_to_product(
  X::PartialFlagVariety, representation::IrrepLevi, offset::Int
)
  weight = _embed_ambient_weight(
    dynkin_type(X), p_dominant_weight(representation), offset
  )
  IrrepLevi(marked_dynkin_type(X), weight)
end

"""
    _lift_bundle_to_product(X, E, offset) -> CompletelyReducibleBundle

Lift a completely reducible bundle from one factor of a product ambient into the
product ambient `X`, embedding each irreducible summand into the coordinate
block starting at `offset`.
"""
function _lift_bundle_to_product(
  X::PartialFlagVariety, E::CompletelyReducibleBundle, offset::Int
)
  _product_factor_range(rank(X), rank(variety(E)), offset)  # validates the fit
  lifted_components = map(components(E)) do component
    _lift_irrep_to_product(X, component, offset)
  end
  CompletelyReducibleBundle(X, lifted_components)
end

# ─── Type-level caches for tangent/cotangent rep lists ───────────────────────
# These avoid repeatedly computing tangent_weights + IrrepLevi decomposition
# for the same MDT.  The caches store only the IrrepLevi component vectors;
# fresh CompletelyReducibleBundle wrappers are created with the correct parent.

const _tangent_reps_cache = let b = _default_cache_budget()
  LRU{MarkedDynkinType,Vector{IrrepLevi}}(;
    maxsize=_cache_maxsize(b, _DEFAULT_STRUCTURAL_FRAC * 0.3),
    by=Base.summarysize,
  )
end
const _cotangent_reps_cache = let b = _default_cache_budget()
  LRU{MarkedDynkinType,Vector{IrrepLevi}}(;
    maxsize=_cache_maxsize(b, _DEFAULT_STRUCTURAL_FRAC * 0.3),
    by=Base.summarysize,
  )
end

function _tangent_reps(mdt::MarkedDynkinType)
  get!(_tangent_reps_cache, mdt) do
    IrrepLevi[IrrepLevi(mdt, w) for w in tangent_weights(mdt)]
  end
end

function _cotangent_reps(mdt::MarkedDynkinType)
  get!(_cotangent_reps_cache, mdt) do
    IrrepLevi[dual(r) for r in _tangent_reps(mdt)]
  end
end

"""
    tangent_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The tangent bundle ``\\mathrm{T}_{G/P}``.

The tangent space at the identity coset decomposes as the direct sum of
root spaces for the positive nonparabolic roots. The tangent bundle
(as an equivariant bundle) is the semisimplification: for each maximal
nonparabolic root class, we get one irreducible Levi summand.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> T = tangent_bundle(X);

julia> n_components(T)
1

julia> rank_bundle(T)
4
```
"""
function tangent_bundle(X::PartialFlagVariety)
  CompletelyReducibleBundle(X, _tangent_reps(marked_dynkin_type(X)))
end

"""
    T(X::PartialFlagVariety) -> CompletelyReducibleBundle

Shorthand for [`tangent_bundle`](@ref): the tangent bundle ``\\mathrm{T}_{G/P}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(T(Gr(2, 4)))
4
```
"""
T(X::PartialFlagVariety) = tangent_bundle(X)

"""
    cotangent_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The cotangent bundle ``\\Omega^1_{G/P} = \\mathrm{T}^\\vee_{G/P}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> rank_bundle(cotangent_bundle(X)) == rank_bundle(tangent_bundle(X))
true
```
"""
function cotangent_bundle(X::PartialFlagVariety)
  CompletelyReducibleBundle(X, _cotangent_reps(marked_dynkin_type(X)))
end

"""
    canonical_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The canonical line bundle ``\\omega_{G/P} = K_{G/P}``.

This is computed directly from the formula

```math
K_{G/P} = -\\sum_{i \\in \\mathrm{marked}} a_i\\,\\omega_i,
\\qquad a_i = \\langle 2(\\rho_G - \\rho_P),\\,\\alpha_i^\\vee\\rangle,
```

without constructing the (co)tangent bundle.  See `anticanonical_degrees`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(1);

julia> K = canonical_bundle(X);

julia> rank_bundle(K)
1

julia> dimensions(K)[1]  # H¹(ℙ¹, 𝒪(-2)) = 1
1

julia> K2 = canonical_bundle(projective_space(4));

julia> dimensions(K2)[4]  # H⁴(ℙ⁴, 𝒪(-5)) = 1
1
```
"""
function canonical_bundle(X::PartialFlagVariety)
  degs = anticanonical_degrees(X)
  line_bundle(X, Vector{Int}(-degs))
end

"""
    anticanonical_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The anticanonical line bundle ``\\omega_{G/P}^\\vee = -K_{G/P}``.

This is computed directly from the formula

```math
-K_{G/P} = \\sum_{i \\in \\mathrm{marked}} a_i\\,\\omega_i,
\\qquad a_i = \\langle 2(\\rho_G - \\rho_P),\\,\\alpha_i^\\vee\\rangle.
```

See `anticanonical_degrees`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(1);

julia> L = anticanonical_bundle(X);

julia> rank_bundle(L)
1

julia> dimensions(L)[0]  # H⁰(ℙ¹, 𝒪(2)) = 3
3

julia> L2 = anticanonical_bundle(projective_space(4));

julia> dimensions(L2)[0]  # H⁰(ℙ⁴, 𝒪(5)) = 126
126
```
"""
function anticanonical_bundle(X::PartialFlagVariety)
  line_bundle(X, anticanonical_degrees(X))
end

"""
    fano_index(X::PartialFlagVariety) -> Int

The Fano index of the partial flag variety ``\\mathrm{G}/\\mathrm{P}``, defined as the gcd of the
anticanonical degrees.

For Picard-rank-1 varieties this is the unique positive integer ``r`` such that

```math
-K_{G/P} = r\\,\\omega_m
```

where ``\\omega_m`` is the ample generator of ``\\mathrm{Pic}(\\mathrm{G}/\\mathrm{P}) \\cong \\mathbb{Z}``.

All partial flag varieties are Fano (the anticanonical bundle is ample), so
the Fano index is always a positive integer.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> fano_index(projective_space(4))
5

julia> fano_index(Gr(2, 5))
5

julia> fano_index(quadric(4))
4
```
"""
function fano_index(X::PartialFlagVariety)
  gcd(anticanonical_degrees(X))
end

_trivial_semisimple_weight(X::PartialFlagVariety) =
  _trivial_semisimple_weight(marked_dynkin_type(X))

# ═══════════════════════════════════════════════════════════════════════════════
#  Monoidal operations on bundles
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dual(E::CompletelyReducibleBundle) -> CompletelyReducibleBundle

The dual bundle ``\\mathcal{E}^\\vee``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> E = tangent_bundle(X);

julia> rank_bundle(dual(E)) == rank_bundle(E)
true
```
"""
function dual(E::CompletelyReducibleBundle)
  CompletelyReducibleBundle(
    variety(E), IrrepLevi[dual(component) for component in E.components]
  )
end

"""
    tensor_product(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)

The tensor product ``\\mathcal{E} \\otimes \\mathcal{F}``.

Uses bilinearity: ``(\\bigoplus_i V_i) \\otimes (\\bigoplus_j \\mathrm{W}_j) = \\bigoplus_{i,j} V_i \\otimes \\mathrm{W}_j``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = tangent_bundle(X);

julia> S = structure_sheaf(X);

julia> rank_bundle(tensor_product(E, S)) == rank_bundle(E)
true
```
"""
function tensor_product(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  X = variety(E)
  Y = variety(F)
  marked_dynkin_type(Y) == marked_dynkin_type(X) || throw(
    ArgumentError(
      "tensor_product requires bundles on the same partial flag variety type."
    ),
  )

  # Deduplicate components to avoid redundant tensor product calls.
  # E.g., O(1)⁶ has 6 identical components → compute ⊗ once, replicate.
  result = IrrepLevi[]
  for (a, mult_a) in _to_counts(E), (b, mult_b) in _to_counts(F)
    product_terms = tensor_product(a, b)
    for _ in 1:(mult_a * mult_b)
      append!(result, product_terms)
    end
  end
  CompletelyReducibleBundle(X, result)
end

"""Multiplicity dict of the irreducible summands of `E`."""
function _to_counts(E::CompletelyReducibleBundle)
  counts = Dict{IrrepLevi,Int}()
  for component in E.components
    counts[component] = get(counts, component, 0) + 1
  end
  counts
end

"""
    exterior_power(E::CompletelyReducibleBundle, k::Int) -> CompletelyReducibleBundle

The k-th exterior power ``\\bigwedge\\nolimits^k \\mathcal{E}``.

For a direct sum ``\\mathcal{E} = \\bigoplus_i V_i``, we use:
``\\bigwedge\\nolimits^k \\mathcal{E} = \\bigoplus_{|\\alpha| = k} \\bigotimes_i \\bigwedge\\nolimits^{\\alpha_i} V_i``

where ``\\alpha`` runs over compositions (multiexponents) of ``k`` in
``n = n_{\\rm components}`` parts.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = tangent_bundle(X);

julia> rank_bundle(exterior_power(E, 0))
1

julia> rank_bundle(exterior_power(E, 1)) == rank_bundle(E)
true
```
"""
function exterior_power(E::CompletelyReducibleBundle, k::Integer)
  k = Int(k)
  r = rank_bundle(E)
  (k < 0 || k > r) && return zero_bundle(E.variety)
  k == 0 && return structure_sheaf(E.variety)
  k == r && return det(E)
  k == 1 && return E

  # Perfect pairing shortcut: ∧^k E ≅ det(E) ⊗ ∧^{r-k} E∨ when k > r/2
  2k > r && return tensor_product(det(E), exterior_power(dual(E), r - k))

  # A group of m copies of an irreducible of rank r absorbs at most m⋅r.
  _power_of_direct_sum(
    exterior_power, E, k; capacity=(c, m) -> Int(fiber_dimension(c)) * m
  )
end

"""
Shared worker for the exterior and symmetric powers of a direct sum: group
equal summands, distribute `k` among the groups, and expand the power of
each ``V^{\\oplus m}`` as

```math
P^{g}(V^{\\oplus m}) = \\bigoplus_{|α| = g} P^{α_1}(V) \\otimes \\cdots \\otimes P^{α_m}(V)
```

for ``P = \\wedge`` or ``\\mathrm{Sym}``, iterating one weakly decreasing
representative per permutation orbit of ``α``, weighted by the multinomial
coefficient counting its permutations.  `capacity(component, mult)` bounds
the exponent a group can absorb (the total rank for exterior powers).
"""
function _power_of_direct_sum(
  power, E::CompletelyReducibleBundle, k::Int; capacity=(c, m) -> typemax(Int)
)
  counts = _to_counts(E)
  unique_comps = collect(keys(counts))
  mults = [counts[component] for component in unique_comps]
  n_groups = length(unique_comps)
  capacities = [capacity(unique_comps[g], mults[g]) for g in 1:n_groups]

  result = IrrepLevi[]

  # Outer multiexponents: how much of k is allocated to each group.
  for group_alloc in multiexponents(n_groups, k)
    any(group_alloc[g] > capacities[g] for g in 1:n_groups) && continue

    group_results = Vector{Vector{Pair{Vector{IrrepLevi},Int}}}()
    skip = false
    for g in 1:n_groups
      component = unique_comps[g]
      group_terms = Pair{Vector{IrrepLevi},Int}[]
      for α in multiexponents(mults[g], group_alloc[g])
        # One canonical (weakly decreasing) representative per orbit.
        issorted(α; rev=true) || continue
        factors = [power(component, a) for a in α]
        any(isempty, factors) && continue
        term = reduce(factors) do acc, factor
          IrrepLevi[t for x in acc for y in factor for t in tensor_product(x, y)]
        end
        push!(group_terms, term => _multinomial_coeff(α))
      end
      if isempty(group_terms)
        skip = true
        break
      end
      push!(group_results, group_terms)
    end
    skip && continue

    # Combine across groups via tensor product.
    _combine_group_results!(result, group_results)
  end

  CompletelyReducibleBundle(E.variety, result)
end

"""
Compute the multinomial coefficient for a multiexponent α:
the number of distinct permutations = n! / (c₁! c₂! ⋯ cₖ!)
where cᵢ counts occurrences of each distinct value.
"""
function _multinomial_coeff(α)
  n = length(α)
  # α is sorted (callers ensure issorted(α; rev=true)); count run-lengths
  result = factorial(big(n))
  i = 1
  @inbounds while i <= n
    j = i
    while j < n && α[j + 1] == α[i]
      j += 1
    end
    c = j - i + 1
    if c > 1
      result = div(result, factorial(big(c)))
    end
    i = j + 1
  end
  Int(result)
end

"""
Combine group results across groups via tensor product, expanding multiplicities.
Each group contributes a list of (IrrepLevi[], multiplicity) pairs.
"""
function _combine_group_results!(result::Vector{IrrepLevi}, group_results)
  if length(group_results) == 1
    for (comps, mult) in group_results[1]
      for _ in 1:mult
        append!(result, comps)
      end
    end
    return nothing
  end

  # Iteratively combine groups
  current_combined = group_results[1]
  for g in 2:length(group_results)
    next_combined = Pair{Vector{IrrepLevi},Int}[]
    for (comps_a, mult_a) in current_combined
      for (comps_b, mult_b) in group_results[g]
        combined = IrrepLevi[]
        for a in comps_a
          for b in comps_b
            append!(combined, tensor_product(a, b))
          end
        end
        push!(next_combined, combined => mult_a * mult_b)
      end
    end
    current_combined = next_combined
  end

  for (comps, mult) in current_combined
    for _ in 1:mult
      append!(result, comps)
    end
  end
end

"""
    symmetric_power(E::CompletelyReducibleBundle, k::Int) -> CompletelyReducibleBundle

The k-th symmetric power ``\\mathrm{Sym}^k \\mathcal{E}``.

For a direct sum ``\\mathcal{E} = \\bigoplus_i V_i``, we use:
``\\mathrm{Sym}^k \\mathcal{E} = \\bigoplus_{|\\alpha| = k} \\bigotimes_i \\mathrm{Sym}^{\\alpha_i} V_i``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = tangent_bundle(X);

julia> rank_bundle(symmetric_power(E, 0))
1

julia> rank_bundle(symmetric_power(E, 1)) == rank_bundle(E)
true
```
"""
function symmetric_power(E::CompletelyReducibleBundle, k::Integer)
  k = Int(k)
  k < 0 && return zero_bundle(E.variety)
  k == 0 && return structure_sheaf(E.variety)
  k == 1 && return E
  _power_of_direct_sum(symmetric_power, E, k)
end

"""
    det(E::CompletelyReducibleBundle) -> CompletelyReducibleBundle
    determinant(E::CompletelyReducibleBundle) -> CompletelyReducibleBundle

The determinant line bundle ``\\det(\\mathcal{E}) = \\bigwedge\\nolimits^{\\mathrm{rk}(\\mathcal{E})} \\mathcal{E}``.

Each summand of `E` is the homogeneous bundle attached to an irreducible
representation ``V`` of the Levi ``L``, and ``\\det`` distributes over direct
sums, ``\\det\\bigl(\\bigoplus_i V_i\\bigr) = \\bigotimes_i \\det V_i``, so the
weights add. For a single ``d``-dimensional irreducible representation ``V`` we
use the factorisation
``L = Z(L)\\cdot [L,L]`` of the reductive Levi into its central torus ``Z(L)``
and semisimple part ``[L,L]``:

- By Schur's lemma ``Z(L)`` acts on ``V`` through a single character ``\\chi``,
  so it acts on ``\\det V = \\bigwedge\\nolimits^d V`` through ``\\chi^d``; the central
  charge of ``\\det V`` is therefore that of ``V`` scaled by ``d``.
- ``[L,L]`` equals its own derived subgroup, so it carries no nontrivial
  characters. The determinant of any ``[L,L]``-representation is therefore
  trivial, and the semisimple highest weight of ``\\det V`` vanishes.

Hence ``\\det V`` is the line bundle with central charge ``d\\,\\chi`` and
trivial semisimple weight, and ``\\det \\mathcal{E}`` is the sum of these over the summands.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> rank_bundle(det(tangent_bundle(X)))
1
```
"""
function det(E::CompletelyReducibleBundle)
  X = variety(E)
  mdt = marked_dynkin_type(X)
  triv_ss = _trivial_semisimple_weight(X)
  λ = WeightLatticeElem(dynkin_type(X))
  for component in components(E)
    fiber_rank = fiber_dimension(component)
    det_rep = IrrepLevi(mdt, fiber_rank .* component.central, triv_ss)
    λ += p_dominant_weight(det_rep)
  end
  CompletelyReducibleBundle(X, λ)
end

"""Alias for [`det`](@ref)."""
determinant(E::CompletelyReducibleBundle) = det(E)

"""
    direct_sum(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)

The direct sum ``\\mathcal{E} \\oplus \\mathcal{F}``.
"""
function direct_sum(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  X = variety(E)
  Y = variety(F)
  marked_dynkin_type(Y) == marked_dynkin_type(X) || throw(
    ArgumentError(
      "direct_sum requires bundles on the same partial flag variety type."
    ),
  )
  CompletelyReducibleBundle(X, vcat(E.components, F.components))
end

Base.:+(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle) = direct_sum(E, F)

# ─── Twist ───────────────────────────────────────────────────────────────────

"""
    twist(E::CompletelyReducibleBundle, i::Int, k::Int=1) -> CompletelyReducibleBundle

Twist ``\\mathcal{E}`` by ``\\mathcal{O}(k)`` at the `i`-th marked node:
``E(k) = \\mathcal{E} \\otimes \\mathcal{L}(k \\omega_{m_i})``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = structure_sheaf(X);

julia> rank_bundle(twist(E, 1, 3))
1
```
"""
function twist(E::CompletelyReducibleBundle, i::Integer, k::Integer=1)
  i, k = Int(i), Int(k)
  mdt = marked_dynkin_type(E.variety)
  marked = marked_nodes(mdt)
  1 <= i <= length(marked) || throw(ArgumentError(
    "Index $i out of range. MDT has $(length(marked)) marked node(s)."
  ))

  m = marked[i]
  ω = fundamental_weight(dynkin_type(mdt), m)
  λ = k * ω
  twist_rep = IrrepLevi(mdt, λ)

  CompletelyReducibleBundle(
    E.variety,
    IrrepLevi[
      t for component in E.components for t in tensor_product(component, twist_rep)
    ],
  )
end

# ─── Arithmetic operators ───────────────────────────────────────────────────

Base.:*(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle) = tensor_product(E, F)

"""
    n * E -> CompletelyReducibleBundle

The `n`-fold direct sum ``\\mathcal{E} \\oplus \\cdots \\oplus \\mathcal{E}`` (n ≥ 1).
Returns the zero bundle when `n == 0` and throws for `n < 0`.
"""
function Base.:*(n::Integer, E::CompletelyReducibleBundle)
  n < 0 && throw(ArgumentError("Cannot multiply a bundle by a negative integer ($n)"))
  n == 0 && return zero_bundle(E.variety)
  CompletelyReducibleBundle(E.variety, repeat(E.components, n))
end
Base.:*(E::CompletelyReducibleBundle, n::Integer) = n * E

const ⊕ = direct_sum
const ⊗ = tensor_product

export ⊕, ⊗

# ─── iszero ──────────────────────────────────────────────────────────

"""
    iszero(E::CompletelyReducibleBundle) -> Bool

Return `true` if `E` is the zero bundle (has no summands).
"""
Base.iszero(E::CompletelyReducibleBundle) = isempty(E.components)

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, E::CompletelyReducibleBundle)
  if isempty(E.components)
    print(io, "0")
  elseif length(E.components) == 1
    print(io, "E", E.components[1])
  else
    parts = ["E" * sprint(show, component) for component in E.components]
    print(io, join(parts, " ⊕ "))
  end
end
