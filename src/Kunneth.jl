# ═══════════════════════════════════════════════════════════════════════════════
#  Structural Künneth recognition
#
#  No product provenance is stored. Recognition therefore runs backwards from
#  product-group representations and ambient presentation terms, accepting only
#  unambiguous external tensor-product decompositions.
# ═══════════════════════════════════════════════════════════════════════════════

# Partition `1:n_vertices` into connected components, joining all vertices in
# one support. Unused vertices remain singleton components.
function _connected_support_components(supports, n_vertices)
  parent = collect(1:n_vertices)
  root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
  for support in supports, vertex in support
    parent[root(vertex)] = root(first(support))
  end
  components = [Int[] for _ in 1:n_vertices]
  for vertex in 1:n_vertices
    push!(components[root(vertex)], vertex)
  end
  filter(!isempty, components)
end

# Split a product-group highest weight at the boundary between the two ambient
# Dynkin types.
function _split_product_irrep(
  irrep::IrrepLevi,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  coefficients = collect(Int, Semisimple.coefficients(p_dominant_weight(irrep)))
  left_type = dynkin_type(left_ambient)
  right_type = dynkin_type(right_ambient)
  left_rank = rank(left_type)
  length(coefficients) == left_rank + rank(right_type) || return nothing

  left_weight = WeightLatticeElem(left_type, coefficients[1:left_rank])
  right_weight = WeightLatticeElem(right_type, coefficients[(left_rank + 1):end])
  (
    IrrepLevi(marked_dynkin_type(left_ambient), left_weight),
    IrrepLevi(marked_dynkin_type(right_ambient), right_weight),
  )
end

# Add `(left, right) => (degree, multiplicity)` to a product grid. Repeated
# factor pairs must have the same degree.
function _add_product_term!(terms, pair, degree::Int, multiplicity::Int=1)
  if haskey(terms, pair)
    old_degree, old_multiplicity = terms[pair]
    old_degree == degree || return false
    terms[pair] = (degree, old_multiplicity + multiplicity)
  else
    terms[pair] = (degree, multiplicity)
  end
  true
end

# Group product-grid keys by shared left or right factor. Returning keys instead
# of sub-dictionaries avoids copying all term values.
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
    [pair for (pair, support) in zip(pairs, supports) if first(support) in block]
    for block in blocks
  ]
end

# Recover factors of one product grid. External products have rectangular
# support, additive degrees, and multiplicative multiplicities.
function _factor_product_terms(terms, component=keys(terms))
  isempty(component) && return nothing
  lefts = unique(first(pair) for pair in component)
  rights = unique(last(pair) for pair in component)
  length(component) == length(lefts) * length(rights) || return nothing

  term_degree(left, right) = first(terms[(left, right)])
  term_multiplicity(left, right) = last(terms[(left, right)])

  # Normalize the opposite degree shift at the first left factor and the common
  # multiplicity scale by the gcd of its row.
  left_anchor, right_anchor = first(lefts), first(rights)
  left_scale = gcd((term_multiplicity(left_anchor, right) for right in rights)...)
  right_data = Dict(
    right => (
      term_degree(left_anchor, right),
      term_multiplicity(left_anchor, right) ÷ left_scale,
    ) for right in rights
  )
  right_scale = last(right_data[right_anchor])
  all(term_multiplicity(left, right_anchor) % right_scale == 0 for left in lefts) ||
    return nothing
  left_data = Dict(
    left => (
      term_degree(left, right_anchor) - term_degree(left_anchor, right_anchor),
      term_multiplicity(left, right_anchor) ÷ right_scale,
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

# Expand factor data into terms grouped by degree. This is shared by filtered
# bundles and zero-locus-bundle presentations.
function _graded_factor_terms(data, ::Type{T}, shift::Int=0) where {T}
  terms = Dict{Int,Vector{T}}()
  for (term, (degree, multiplicity)) in data
    append!(get!(terms, degree + shift, T[]), fill(term, multiplicity))
  end
  terms
end

# Rebuild a filtered factor from its irreducible factor data.
function _filtered_bundle_from_factor_data(ambient::PartialFlagVariety, data)
  terms = _graded_factor_terms(data, IrrepLevi)
  FilteredBundle(
    ambient,
    CompletelyReducibleBundle[
      CompletelyReducibleBundle(ambient, terms[degree]) for
      degree in sort!(collect(keys(terms)))
    ],
  )
end

# Split every irreducible summand of a completely reducible ambient bundle.
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

# Normalize a one-step filtration to its total bundle. For a genuine filtration,
# factor its irreducible grid and verify the reconstructed external product.
function _factor_ambient_bundle_on_product(
  bundle::FilteredBundle,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  n_filtration_steps(bundle) == 1 && return _factor_ambient_bundle_on_product(
    total_bundle(bundle), left_ambient, right_ambient
  )

  terms = Dict{Tuple{IrrepLevi,IrrepLevi},Tuple{Int,Int}}()
  for (degree, piece) in enumerate(graded_pieces(bundle)), irrep in components(piece)
    split = _split_product_irrep(irrep, left_ambient, right_ambient)
    split === nothing && return nothing
    _add_product_term!(terms, split, degree) || return nothing
  end
  factor_data = _factor_product_terms(terms)
  factor_data === nothing && return nothing

  left_data, right_data = factor_data
  left = _filtered_bundle_from_factor_data(left_ambient, left_data)
  right = _filtered_bundle_from_factor_data(right_ambient, right_data)
  normalize_factor(F) = n_filtration_steps(F) == 1 ? total_bundle(F) : F
  factor_pair = normalize_factor(left), normalize_factor(right)
  external_tensor_product(factor_pair...) == bundle ? [factor_pair] : nothing
end

# Rebuild a zero-locus-bundle factor with the chosen cohomological shift.
function _presentation_from_factor_data(
  locus::ZeroLocus, data, degree_shift::Int
)
  terms = _graded_factor_terms(data, _AmbientBundle, degree_shift)
  ZeroLocusBundle(locus, _AmbientBundlePresentation(terms))
end

# Factor one connected presentation component and choose the unique shift for
# which both factor presentations have positive rank in degree zero.
function _factor_presentation_component(
  terms, component, left_locus::ZeroLocus, right_locus::ZeroLocus
)
  factor_data = _factor_product_terms(terms, component)
  factor_data === nothing && return nothing
  left_data, right_data = factor_data

  shifts = intersect(
    Set(-degree for (_, (degree, _)) in left_data),
    Set(degree for (_, (degree, _)) in right_data),
  )
  factor_rank(data, shift) = sum(
    (isodd(degree + shift) ? -1 : 1) * multiplicity * rank(bundle) for
    (bundle, (degree, multiplicity)) in data;
    init=0,
  )
  filter!(
    shift -> factor_rank(left_data, shift) > 0 && factor_rank(right_data, -shift) > 0,
    shifts,
  )
  length(shifts) == 1 || return nothing

  shift = only(shifts)
  (
    _presentation_from_factor_data(left_locus, left_data, shift),
    _presentation_from_factor_data(right_locus, right_data, -shift),
  )
end

# Recognize a presentation as a direct sum of external products for one fixed
# bipartition of its zero locus.
function _kunneth_decomposition(
  F::ZeroLocusBundle, left_locus::ZeroLocus, right_locus::ZeroLocus
)
  left_ambient = ambient_variety(left_locus)
  right_ambient = ambient_variety(right_locus)
  terms = Dict{Tuple{_AmbientBundle,_AmbientBundle},Tuple{Int,Int}}()

  for (degree, summands) in F.presentation.terms, summand in summands
    factorizations = _factor_ambient_bundle_on_product(summand, left_ambient, right_ambient)
    factorizations === nothing && return nothing
    for (left, right) in factorizations
      _add_product_term!(terms, (left, right), degree) || return nothing
    end
  end
  isempty(terms) && return nothing

  decomposition = Tuple{ZeroLocusBundle,ZeroLocusBundle}[]
  for component in _product_term_components(terms)
    factor_pair = _factor_presentation_component(terms, component, left_locus, right_locus)
    factor_pair === nothing && return nothing
    push!(decomposition, factor_pair)
  end
  decomposition
end

# Try each contiguous bipartition of the recognized zero-locus factors.
function _kunneth_decomposition(F::ZeroLocusBundle)
  Z = variety(F)
  locus_factors = factors(Z)
  for split in 1:(length(locus_factors) - 1)
    left = reduce(product, locus_factors[1:split])
    right = reduce(product, locus_factors[(split + 1):end])
    product(left, right) == Z || continue
    decomposition = _kunneth_decomposition(F, left, right)
    decomposition === nothing || return decomposition
  end
  nothing
end

# Convolve determined factor cohomologies; otherwise leave the generic
# presentation solver to retain its symbolic answer.
function _kunneth_cohomology(F::ZeroLocusBundle)
  decomposition = _kunneth_decomposition(F)
  decomposition === nothing && return nothing

  d = dimension(variety(F))
  entries = zeros(BigInt, d + 1)
  for (left, right) in decomposition
    left_cohomology = cohomology(left)
    right_cohomology = cohomology(right)
    is_determined(left_cohomology) && is_determined(right_cohomology) || return nothing
    for p in 0:left_cohomology.max_degree, q in 0:right_cohomology.max_degree
      entries[p + q + 1] +=
        left_cohomology[p].constant * right_cohomology[q].constant
    end
  end
  Cohomology{AffineExpr}(AffineExpr.(entries), d)
end
