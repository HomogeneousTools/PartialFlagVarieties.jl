# ═══════════════════════════════════════════════════════════════════════════════
#  External products of ambient bundles
# ═══════════════════════════════════════════════════════════════════════════════

export external_direct_sum, external_tensor_product

const _AmbientBundle = Union{CompletelyReducibleBundle,FilteredBundle}

# Lift every graded piece of `F` into one factor of `product_ambient`.
function _lift_bundle_to_product(
  product_ambient::PartialFlagVariety, F::FilteredBundle, offset::Int
)
  FilteredBundle(
    product_ambient,
    [
      _lift_bundle_to_product(product_ambient, piece, offset) for
      piece in graded_pieces(F)
    ],
  )
end

# Lift two ambient bundles to the common product ambient, placing the right
# bundle after the left Dynkin block.
function _lift_external_factors(
  product_ambient::PartialFlagVariety,
  E::_AmbientBundle,
  F::_AmbientBundle,
  right_offset::Int,
)
  (
    _lift_bundle_to_product(product_ambient, E, 0),
    _lift_bundle_to_product(product_ambient, F, right_offset),
  )
end

"""
    external_tensor_product(E, F) -> Bundle

Return the external tensor product
``p_X^*\\mathcal{E} \\otimes p_Y^*\\mathcal{F}`` on ``X \\times Y``, where
`E` lives on `X` and `F` lives on `Y`.

The arguments may be [`CompletelyReducibleBundle`](@ref)s or
[`FilteredBundle`](@ref)s. If either argument is filtered, the result carries
the induced tensor-product filtration.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X, Y = projective_space(1), projective_space(2);

julia> L = external_tensor_product(line_bundle(X, 2), line_bundle(Y, 3));

julia> L == line_bundle(product(X, Y), [2, 3])
true
```
"""
function external_tensor_product(E::_AmbientBundle, F::_AmbientBundle)
  product_ambient = product(variety(E), variety(F))
  right_offset = rank(dynkin_type(variety(E)))
  lifted_E, lifted_F = _lift_external_factors(product_ambient, E, F, right_offset)
  tensor_product(lifted_E, lifted_F)
end

"""
    external_direct_sum(
      E::CompletelyReducibleBundle,
      F::CompletelyReducibleBundle,
    ) -> CompletelyReducibleBundle

Return the external direct sum
``p_X^*\\mathcal{E} \\oplus p_Y^*\\mathcal{F}`` on ``X \\times Y``, where
`E` lives on `X` and `F` lives on `Y`.

This method is restricted to completely reducible bundles. A
[`FilteredBundle`](@ref) records only an ordered list of graded pieces, without
filtration indices that could canonically align the two summand filtrations.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X, Y = projective_space(1), projective_space(2);

julia> E = external_direct_sum(line_bundle(X, 2), line_bundle(Y, 3));

julia> E == direct_sum(line_bundle(product(X, Y), [2, 0]),
                       line_bundle(product(X, Y), [0, 3]))
true
```
"""
function external_direct_sum(
  E::CompletelyReducibleBundle, F::CompletelyReducibleBundle
)
  product_ambient = product(variety(E), variety(F))
  right_offset = rank(dynkin_type(variety(E)))
  lifted_E, lifted_F = _lift_external_factors(product_ambient, E, F, right_offset)
  direct_sum(lifted_E, lifted_F)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Structural factorization of ambient external products
# ═══════════════════════════════════════════════════════════════════════════════

# Partition `1:n` into connected components, joining all vertices occurring in
# one support. Unused vertices remain singleton components.
function _connected_support_components(supports, n)
  parent = collect(1:n)
  root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
  for support in supports, vertex in support
    parent[root(vertex)] = root(first(support))
  end
  components = [Int[] for _ in 1:n]
  for vertex in 1:n
    push!(components[root(vertex)], vertex)
  end
  filter(!isempty, components)
end

# Split the highest weight of a product-group representation into its two
# Dynkin blocks. This depends only on the requested ambient bipartition.
function _split_product_irrep(
  irrep::IrrepLevi,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  coefficients = collect(Int, Semisimple.coefficients(p_dominant_weight(irrep)))
  left_type = dynkin_type(left_ambient)
  right_type = dynkin_type(right_ambient)
  left_rank = rank(left_type)
  right_rank = rank(right_type)
  length(coefficients) == left_rank + right_rank || return nothing

  left_weight = WeightLatticeElem(left_type, coefficients[1:left_rank])
  right_weight = WeightLatticeElem(right_type, coefficients[(left_rank + 1):end])
  (
    IrrepLevi(marked_dynkin_type(left_ambient), left_weight),
    IrrepLevi(marked_dynkin_type(right_ambient), right_weight),
  )
end

# Product terms are stored as `(left, right) => (degree, multiplicity)`.
# Repeated pairs must occur in the same degree to admit this factorization.
function _add_product_term!(terms, left, right, degree::Int, multiplicity::Int=1)
  key = (left, right)
  if haskey(terms, key)
    old_degree, old_multiplicity = terms[key]
    old_degree == degree || return false
    terms[key] = (degree, old_multiplicity + multiplicity)
  else
    terms[key] = (degree, multiplicity)
  end
  true
end

# Partition pairs into candidate summands. Pairs are linked when they share a
# left or right factor; the transitive closure gives the maximal groups.
function _product_term_components(terms)
  pairs = collect(keys(terms))
  lefts = unique(first(pair) for pair in pairs)
  rights = unique(last(pair) for pair in pairs)
  left_index = Dict(left => index for (index, left) in enumerate(lefts))
  right_index = Dict(right => index for (index, right) in enumerate(rights))
  supports = [
    (left_index[left], length(lefts) + right_index[right]) for
    (left, right) in pairs
  ]
  blocks = _connected_support_components(supports, length(lefts) + length(rights))
  [
    typeof(terms)(
      pair => terms[pair] for (pair, support) in zip(pairs, supports) if
      first(support) in block
    ) for block in blocks
  ]
end

# Recover two factors from a connected term group. External products have
# rectangular support, additive degrees, and multiplicative multiplicities.
function _factor_product_terms(terms)
  isempty(terms) && return nothing
  lefts = unique(first(pair) for pair in keys(terms))
  rights = unique(last(pair) for pair in keys(terms))
  length(terms) == length(lefts) * length(rights) || return nothing

  degree(left, right) = first(terms[(left, right)])
  multiplicity(left, right) = last(terms[(left, right)])

  # The split is defined up to an opposite degree shift and a common integer
  # factor. Normalize with left degree zero and a primitive right first row.
  left0, right0 = first(lefts), first(rights)
  left_scale = gcd((multiplicity(left0, right) for right in rights)...)
  right_data = Dict(
    right => (
      degree(left0, right),
      multiplicity(left0, right) ÷ left_scale,
    ) for right in rights
  )
  right_scale = last(right_data[right0])
  all(multiplicity(left, right0) % right_scale == 0 for left in lefts) ||
    return nothing
  left_data = Dict(
    left => (
      degree(left, right0) - degree(left0, right0),
      multiplicity(left, right0) ÷ right_scale,
    ) for left in lefts
  )

  for left in lefts, right in rights
    left_degree, left_multiplicity = left_data[left]
    right_degree, right_multiplicity = right_data[right]
    terms[(left, right)] == (
      left_degree + right_degree,
      left_multiplicity * right_multiplicity,
    ) || return nothing
  end
  left_data, right_data
end

# Rebuild a filtered factor from `(filtration degree, multiplicity)` data for
# each irreducible representation.
function _filtered_bundle_from_factor_data(ambient::PartialFlagVariety, data)
  pieces = Dict{Int,Vector{IrrepLevi}}()
  for (irrep, (degree, multiplicity)) in data
    append!(get!(pieces, degree, IrrepLevi[]), fill(irrep, multiplicity))
  end
  ordered_pieces = CompletelyReducibleBundle[
    CompletelyReducibleBundle(ambient, pieces[degree]) for
    degree in sort!(collect(keys(pieces)))
  ]
  FilteredBundle(ambient, ordered_pieces)
end

# Factor a multi-step filtered bundle by factoring its grid of irreducible
# graded pieces, then verify the reconstructed external product structurally.
function _factor_filtered_bundle_on_product(
  bundle::FilteredBundle,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  terms = Dict{Tuple{IrrepLevi,IrrepLevi},Tuple{Int,Int}}()
  for (degree, piece) in enumerate(graded_pieces(bundle))
    for irrep in components(piece)
      split = _split_product_irrep(irrep, left_ambient, right_ambient)
      split === nothing && return nothing
      _add_product_term!(terms, split[1], split[2], degree) || return nothing
    end
  end
  factor_data = _factor_product_terms(terms)
  factor_data === nothing && return nothing
  left_data, right_data = factor_data

  left = _filtered_bundle_from_factor_data(left_ambient, left_data)
  right = _filtered_bundle_from_factor_data(right_ambient, right_data)
  normalize(F) = n_filtration_steps(F) == 1 ? total_bundle(F) : F
  factor_pair = normalize(left), normalize(right)
  external_tensor_product(factor_pair...) == bundle ? factor_pair : nothing
end

# Split every irreducible summand of a completely reducible ambient bundle.
# The result is its direct-sum decomposition into external tensor products.
function _factor_ambient_bundle_on_product(
  bundle::CompletelyReducibleBundle,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  splits = [
    _split_product_irrep(irrep, left_ambient, right_ambient) for
    irrep in components(bundle)
  ]
  any(isnothing, splits) && return nothing
  [
    (
      CompletelyReducibleBundle(left_ambient, IrrepLevi[split[1]]),
      CompletelyReducibleBundle(right_ambient, IrrepLevi[split[2]]),
    ) for split in splits
  ]
end

# A one-step filtration is normalized to its total bundle; a genuine
# multi-step filtration is factored as one filtered external product.
function _factor_ambient_bundle_on_product(
  bundle::FilteredBundle,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  n_filtration_steps(bundle) == 1 && return _factor_ambient_bundle_on_product(
    total_bundle(bundle), left_ambient, right_ambient
  )
  factorization = _factor_filtered_bundle_on_product(bundle, left_ambient, right_ambient)
  factorization === nothing ? nothing : [factorization]
end
