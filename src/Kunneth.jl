# ═══════════════════════════════════════════════════════════════════════════════
#  Structural Künneth recognition
#
#  No product provenance is stored. Recognition therefore runs backwards from
#  product-group representations and ambient presentation terms, accepting only
#  unambiguous external tensor-product decompositions.
# ═══════════════════════════════════════════════════════════════════════════════

# Split a product-group highest weight at the boundary between the two ambient
# Dynkin types. Callers have already validated the product ambient.
function _split_product_irrep(
  irrep::IrrepLevi,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  coefficients = collect(Int, Semisimple.coefficients(p_dominant_weight(irrep)))
  left_rank = rank(dynkin_type(left_ambient))
  (
    IrrepLevi(marked_dynkin_type(left_ambient), coefficients[1:left_rank]),
    IrrepLevi(marked_dynkin_type(right_ambient), coefficients[(left_rank + 1):end]),
  )
end

# Add `(left, right) => (degree, multiplicity)` to a product grid. Repeated
# factor pairs must have the same degree.
function _add_product_term!(terms, pair, degree::Int, multiplicity::Int=1)
  old_degree, old_multiplicity = get(terms, pair, (degree, 0))
  old_degree == degree || return false
  terms[pair] = (degree, old_multiplicity + multiplicity)
  true
end

# Group product-grid keys by shared left or right factor. Presentations are
# small, so a direct component walk is clearer than a separate graph structure.
function _product_term_components(terms)
  remaining = Set(keys(terms))
  components = Vector{Vector{keytype(terms)}}()
  while !isempty(remaining)
    component = keytype(terms)[pop!(remaining)]
    cursor = 1
    while cursor <= length(component)
      pair = component[cursor]
      neighbors = [
        candidate for candidate in remaining if
        first(candidate) == first(pair) || last(candidate) == last(pair)
      ]
      append!(component, neighbors)
      setdiff!(remaining, neighbors)
      cursor += 1
    end
    push!(components, component)
  end
  components
end

# Recover factors of one product grid. External products have rectangular
# support, additive degrees, and multiplicative multiplicities.
function _factor_product_terms(terms, component=keys(terms))
  isempty(component) && return nothing
  lefts = unique(first(pair) for pair in component)
  rights = unique(last(pair) for pair in component)
  length(component) == length(lefts) * length(rights) || return nothing

  degrees = [first(terms[(left, right)]) for left in lefts, right in rights]
  multiplicities = [last(terms[(left, right)]) for left in lefts, right in rights]
  degrees == degrees[:, 1] .+ degrees[1, :]' .- degrees[1, 1] || return nothing
  multiplicities .* multiplicities[1, 1] ==
  multiplicities[:, 1] * multiplicities[1, :]' || return nothing

  # Normalize the opposite degree shift at the first left factor and the common
  # multiplicity scale by the gcd of its row. The rank-one identity makes the
  # division of the first column below exact.
  row_gcd = reduce(gcd, @view multiplicities[1, :])
  right_multiplicities = multiplicities[1, :] .÷ row_gcd
  left_multiplicities = multiplicities[:, 1] .÷ first(right_multiplicities)
  left_data = Dict(
    zip(lefts, zip(degrees[:, 1] .- degrees[1, 1], left_multiplicities))
  )
  right_data = Dict(zip(rights, zip(degrees[1, :], right_multiplicities)))
  left_data, right_data
end

# Expand factor data into terms grouped by degree. This is shared by filtered
# bundles and zero-locus-bundle presentations.
function _graded_factor_terms(data::AbstractDict{T}, shift::Int=0) where {T}
  terms = Dict{Int,Vector{T}}()
  for (term, (degree, multiplicity)) in data
    append!(get!(terms, degree + shift, T[]), fill(term, multiplicity))
  end
  terms
end

# Split every irreducible summand of a completely reducible ambient bundle.
function _factor_ambient_bundle_on_product(
  bundle::CompletelyReducibleBundle,
  left_ambient::PartialFlagVariety,
  right_ambient::PartialFlagVariety,
)
  map(components(bundle)) do irrep
    left, right = _split_product_irrep(irrep, left_ambient, right_ambient)
    (
      CompletelyReducibleBundle(left_ambient, [left]),
      CompletelyReducibleBundle(right_ambient, [right]),
    )
  end
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
    _add_product_term!(terms, split, degree) || return nothing
  end
  factor_data = _factor_product_terms(terms)
  factor_data === nothing && return nothing

  factor_pair = map((left_ambient, right_ambient), factor_data) do ambient, data
    factor = _filtered_bundle_from_graded_components(
      ambient, _graded_factor_terms(data)
    )
    n_filtration_steps(factor) == 1 ? total_bundle(factor) : factor
  end
  external_tensor_product(factor_pair...) == bundle ? [factor_pair] : nothing
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
  map((left_locus, right_locus), factor_data, (shift, -shift)) do locus, data, offset
    ZeroLocusBundle(
      locus, _AmbientBundlePresentation(_graded_factor_terms(data, offset))
    )
  end
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
