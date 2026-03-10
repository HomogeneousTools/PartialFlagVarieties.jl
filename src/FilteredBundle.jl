# ═══════════════════════════════════════════════════════════════════════════════
#  FilteredBundle — equivariant bundle with a filtration on G/P
#
#  A filtered bundle is an equivariant bundle equipped with a filtration
#  by equivariant subbundles. It is encoded by its associated graded pieces,
#  each of which is a CompletelyReducibleBundle. The ordering of the pieces
#  records the filtration.
#
#  The main example is the tangent bundle T_{G/P}, which has a natural
#  filtration by root height: the associated graded pieces are indexed
#  by the height levels of the positive nonparabolic roots, and each
#  piece is the sum of the Levi representations at that height.
#  (cf. Lemma 2.1 of arXiv:1606.04076)
# ═══════════════════════════════════════════════════════════════════════════════

export FilteredBundle
export graded_pieces, total_bundle, filtered_tangent_bundle, n_filtration_steps

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    FilteredBundle{MDT}

An equivariant vector bundle on ``G/P`` equipped with a filtration by
equivariant subbundles. Stored as an ordered list of associated graded
pieces (each a [`CompletelyReducibleBundle`](@ref)).

The filtration is encoded by the ordering: `graded_pieces(F)[1]` is the
bottom piece (smallest filtration step), and `graded_pieces(F)[end]` is
the top piece.

# Fields
- `variety::PartialFlagVariety{MDT}`: the partial flag variety
- `pieces::Vector{CompletelyReducibleBundle{MDT}}`: the associated graded pieces

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> n_filtration_steps(F)
1
```
"""
struct FilteredBundle{MDT<:MarkedDynkinType} <: Bundle{MDT}
  variety::PartialFlagVariety{MDT}
  pieces::Vector{CompletelyReducibleBundle{MDT}}

  function FilteredBundle{MDT}(
    variety::PartialFlagVariety{MDT},
    pieces::Vector{CompletelyReducibleBundle{MDT}}
  ) where {MDT}
    new{MDT}(variety, pieces)
  end
end

function FilteredBundle(
  variety::PartialFlagVariety{MDT},
  pieces::Vector{CompletelyReducibleBundle{MDT}}
) where {MDT}
  FilteredBundle{MDT}(variety, pieces)
end

# ─── Accessors ───────────────────────────────────────────────────────────────

"""
    variety(F::FilteredBundle) -> PartialFlagVariety

Return the partial flag variety on which this filtered bundle lives.
"""
variety(F::FilteredBundle) = F.variety

"""
    graded_pieces(F::FilteredBundle) -> Vector{CompletelyReducibleBundle}

Return the associated graded pieces of the filtration, ordered from
bottom to top.
"""
graded_pieces(F::FilteredBundle) = F.pieces

"""
    n_filtration_steps(F::FilteredBundle) -> Int

Return the number of filtration steps (graded pieces).
"""
n_filtration_steps(F::FilteredBundle) = length(F.pieces)

"""
    total_bundle(F::FilteredBundle) -> CompletelyReducibleBundle

Return the total bundle as a [`CompletelyReducibleBundle`](@ref),
forgetting the filtration (i.e., the direct sum of all graded pieces).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> rank_bundle(total_bundle(F)) == dimension(X)
true
```
"""
function total_bundle(F::FilteredBundle{MDT}) where {MDT}
  all_comps = IrrepLevi{MDT}[]
  for piece in F.pieces
    append!(all_comps, components(piece))
  end
  CompletelyReducibleBundle{MDT}(F.variety, all_comps)
end

"""
    rank_bundle(F::FilteredBundle) -> BigInt

Return the total rank of the filtered bundle.
"""
function rank_bundle(F::FilteredBundle)
  sum(rank_bundle(p) for p in F.pieces; init=BigInt(0))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Tensor product of FilteredBundle with CompletelyReducibleBundle
# ═══════════════════════════════════════════════════════════════════════════════

"""
    tensor_product(F::FilteredBundle, E::CompletelyReducibleBundle) -> FilteredBundle

Tensor each graded piece of `F` with `E`, preserving the filtration.
"""
function tensor_product(F::FilteredBundle, E::CompletelyReducibleBundle)
  new_pieces = [tensor_product(p, E) for p in F.pieces]
  FilteredBundle(F.variety, new_pieces)
end

function tensor_product(E::CompletelyReducibleBundle, F::FilteredBundle)
  tensor_product(F, E)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Filtered tangent bundle
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _root_heights(::Type{MDT}) -> Vector{Tuple{Int, Vector{WeightLatticeElem}}}

Compute the nonparabolic root height decomposition. Returns pairs
`(height, weights)` sorted by increasing height, where `weights` are
the highest weights of the Levi representations appearing at that height
level.

The nonparabolic height of a positive root ``\\alpha = \\sum_i a_i \\alpha_i``
is ``\\sum_{j \\in I} a_j`` where ``I`` is the set of marked (crossed-out)
nodes.
"""
function _root_heights(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  R = rank(DT)
  RS = RootSystem(DT)
  Marked = marked_nodes(MDT)
  unmarked = unmarked_nodes(MDT)

  nonpar_roots = positive_nonparabolic_roots(MDT)

  # Group roots by nonparabolic height
  height_groups = Dict{Int,Vector{typeof(first(nonpar_roots))}}()
  for α in nonpar_roots
    h = _nonparabolic_height(coefficients(α), Marked)
    if !haskey(height_groups, h)
      height_groups[h] = typeof(α)[]
    end
    push!(height_groups[h], α)
  end

  # Simple parabolic roots (at unmarked positions)
  simple_par = [simple_root(RS, i) for i in unmarked]

  # For each height level, find maximal roots (= highest weights for the Levi action)
  result = Tuple{Int,Vector{WeightLatticeElem}}[]

  for h in sort(collect(keys(height_groups)))
    roots_at_h = height_groups[h]
    roots_set = Set(coefficients(α) for α in roots_at_h)

    # A root α at height h is maximal if α + β is not a root at height h
    # for any simple parabolic root β
    maximal = filter(roots_at_h) do α
      !any(s -> (coefficients(α) + coefficients(s)) in roots_set, simple_par)
    end

    push!(result, (h, [WeightLatticeElem(α) for α in maximal]))
  end

  result
end

"""
    filtered_tangent_bundle(X::PartialFlagVariety) -> FilteredBundle

The tangent bundle ``T_{G/P}`` with its natural filtration by
nonparabolic root height.

The associated graded pieces are indexed by the root height levels
``h = 1, 2, \\ldots``, and each piece contains the irreducible Levi
representations whose highest weight has nonparabolic height ``h``.

This is the filtration described in Lemma 2.1 of
[arXiv:1606.04076](https://arxiv.org/abs/1606.04076).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> n_filtration_steps(F)
1

julia> rank_bundle(F) == dimension(X)
true
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(2);

julia> F = filtered_tangent_bundle(X);

julia> n_filtration_steps(F)
1
```
"""
function filtered_tangent_bundle(X::PartialFlagVariety{MDT}) where {MDT<:MarkedDynkinType}
  height_data = _root_heights(MDT)

  pieces = CompletelyReducibleBundle{MDT}[]
  for (h, ws) in height_data
    reps = [IrrepLevi(MDT, w) for w in ws]
    push!(pieces, CompletelyReducibleBundle{MDT}(X, reps))
  end

  FilteredBundle(X, pieces)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Exterior and symmetric powers of a FilteredBundle
#
#  Given F with graded pieces gr_1, …, gr_s, the k-th exterior power ∧^k F
#  has an induced filtration whose associated graded is
#
#      ⊕_{|α|=k} ∧^{α_1} gr_1 ⊗ ⋯ ⊗ ∧^{α_s} gr_s
#
#  where the sum runs over all multiexponents α = (α_1,…,α_s) with |α| = k.
#  Each such term has "filtration weight" Σ_i i⋅α_i (i.e., the position in
#  the filtration is determined by how deep in the filtration the factors come
#  from). Terms with the same filtration weight form a single graded piece.
#
#  This is the standard construction; cf. Lemma 23 of arXiv:1911.09414 for
#  the specific application to the tangent bundle's height filtration on G/P.
#
#  Symmetric powers work analogously with Sym replacing ∧.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    exterior_power(F::FilteredBundle, k::Integer) -> FilteredBundle

The ``k``-th exterior power ``\\bigwedge^k F`` of a filtered bundle, equipped
with the induced filtration.

The associated graded of ``\\bigwedge^k F`` is
```math
\\bigoplus_{|\\alpha|=k} \\bigwedge^{\\alpha_1} \\mathrm{gr}_1 \\otimes \\cdots \\otimes \\bigwedge^{\\alpha_s} \\mathrm{gr}_s
```
where ``\\mathrm{gr}_i`` are the graded pieces of ``F``, and the filtration
is ordered by filtration weight ``\\sum_i i \\cdot \\alpha_i``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> rank_bundle(exterior_power(F, 2)) == binomial(dimension(X), 2)
true
```
"""
function exterior_power(F::FilteredBundle{MDT}, k::Integer) where {MDT}
  k = Int(k)
  s = n_filtration_steps(F)
  k < 0 && return FilteredBundle(F.variety, CompletelyReducibleBundle{MDT}[])
  k == 0 && return FilteredBundle(F.variety, [structure_sheaf(F.variety)])
  k == 1 && return F

  pieces = graded_pieces(F)
  ranks = [Int(rank_bundle(p)) for p in pieces]

  # Collect terms by filtration weight
  weight_terms = Dict{Int,Vector{CompletelyReducibleBundle{MDT}}}()

  for α in multiexponents(s, k)
    # Skip if any α_i exceeds the rank of gr_i
    any(α[i] > ranks[i] for i in 1:s) && continue

    # Compute ∧^{α_1} gr_1 ⊗ ⋯ ⊗ ∧^{α_s} gr_s
    factors = CompletelyReducibleBundle{MDT}[]
    skip = false
    for i in 1:s
      w_i = exterior_power(pieces[i], α[i])
      if iszero(w_i)
        skip = true
        break
      end
      push!(factors, w_i)
    end
    skip && continue

    # Tensor all factors together
    term = factors[1]
    for i in 2:length(factors)
      term = tensor_product(term, factors[i])
    end

    # Filtration weight = Σ i * α_i (1-indexed)
    fw = sum(i * α[i] for i in 1:s)
    if !haskey(weight_terms, fw)
      weight_terms[fw] = CompletelyReducibleBundle{MDT}[]
    end
    push!(weight_terms[fw], term)
  end

  # Assemble graded pieces ordered by filtration weight
  result_pieces = CompletelyReducibleBundle{MDT}[]
  for fw in sort(collect(keys(weight_terms)))
    # Direct sum of all terms at this filtration weight
    all_comps = IrrepLevi{MDT}[]
    for t in weight_terms[fw]
      append!(all_comps, components(t))
    end
    push!(result_pieces, CompletelyReducibleBundle{MDT}(F.variety, all_comps))
  end

  FilteredBundle(F.variety, result_pieces)
end

"""
    symmetric_power(F::FilteredBundle, k::Integer) -> FilteredBundle

The ``k``-th symmetric power ``\\mathrm{Sym}^k F`` of a filtered bundle,
equipped with the induced filtration.

The associated graded of ``\\mathrm{Sym}^k F`` is
```math
\\bigoplus_{|\\alpha|=k} \\mathrm{Sym}^{\\alpha_1} \\mathrm{gr}_1 \\otimes \\cdots \\otimes \\mathrm{Sym}^{\\alpha_s} \\mathrm{gr}_s
```

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> rank_bundle(symmetric_power(F, 2)) == binomial(dimension(X) + 1, 2)
true
```
"""
function symmetric_power(F::FilteredBundle{MDT}, k::Integer) where {MDT}
  k = Int(k)
  s = n_filtration_steps(F)
  k < 0 && return FilteredBundle(F.variety, CompletelyReducibleBundle{MDT}[])
  k == 0 && return FilteredBundle(F.variety, [structure_sheaf(F.variety)])
  k == 1 && return F

  pieces = graded_pieces(F)

  # Collect terms by filtration weight
  weight_terms = Dict{Int,Vector{CompletelyReducibleBundle{MDT}}}()

  for α in multiexponents(s, k)
    # Compute Sym^{α_1} gr_1 ⊗ ⋯ ⊗ Sym^{α_s} gr_s
    factors = CompletelyReducibleBundle{MDT}[]
    skip = false
    for i in 1:s
      s_i = symmetric_power(pieces[i], α[i])
      if iszero(s_i)
        skip = true
        break
      end
      push!(factors, s_i)
    end
    skip && continue

    term = factors[1]
    for i in 2:length(factors)
      term = tensor_product(term, factors[i])
    end

    fw = sum(i * α[i] for i in 1:s)
    if !haskey(weight_terms, fw)
      weight_terms[fw] = CompletelyReducibleBundle{MDT}[]
    end
    push!(weight_terms[fw], term)
  end

  result_pieces = CompletelyReducibleBundle{MDT}[]
  for fw in sort(collect(keys(weight_terms)))
    all_comps = IrrepLevi{MDT}[]
    for t in weight_terms[fw]
      append!(all_comps, components(t))
    end
    push!(result_pieces, CompletelyReducibleBundle{MDT}(F.variety, all_comps))
  end

  FilteredBundle(F.variety, result_pieces)
end

"""
    dual(F::FilteredBundle) -> FilteredBundle

The dual of a filtered bundle. The filtration is reversed: if ``F`` has
pieces ``\\mathrm{gr}_1, \\ldots, \\mathrm{gr}_s`` (bottom to top), then
``F^\\vee`` has pieces ``\\mathrm{gr}_s^\\vee, \\ldots, \\mathrm{gr}_1^\\vee``
(the dual reverses the filtration order).
"""
function dual(F::FilteredBundle{MDT}) where {MDT}
  FilteredBundle(F.variety, [dual(p) for p in reverse(F.pieces)])
end

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, F::FilteredBundle{MDT}) where {MDT}
  n = n_filtration_steps(F)
  r = rank_bundle(F)
  print(io, "FilteredBundle(rank $r, $n layer(s))")
end
