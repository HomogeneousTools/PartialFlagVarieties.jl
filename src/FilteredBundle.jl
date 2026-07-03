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
  FilteredBundle

An equivariant vector bundle on ``G/P`` equipped with a filtration by
equivariant subbundles. Stored as an ordered list of associated graded
pieces (each a [`CompletelyReducibleBundle`](@ref)).

The filtration is encoded by the ordering: `graded_pieces(F)[1]` is the
bottom piece (smallest filtration step), and `graded_pieces(F)[end]` is
the top piece.

# Fields
- `variety::PartialFlagVariety`: the partial flag variety
- `pieces::Vector{CompletelyReducibleBundle}`: the associated graded pieces

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> n_filtration_steps(F)
1
```
"""
struct FilteredBundle <: Bundle
  variety::PartialFlagVariety
  pieces::Vector{CompletelyReducibleBundle}

  function FilteredBundle(
    X::PartialFlagVariety,
    pieces::Vector{CompletelyReducibleBundle},
  )
    for (idx, piece) in enumerate(pieces)
      marked_dynkin_type(variety(piece)) == marked_dynkin_type(X) || throw(
        ArgumentError(
          "Filtered bundle piece $idx lives on $(variety(piece)), expected $X."
        ),
      )
    end
    new(X, pieces)
  end
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
function total_bundle(F::FilteredBundle)
  all_components = IrrepLevi[c for piece in F.pieces for c in components(piece)]
  CompletelyReducibleBundle(F.variety, all_components)
end

"""
    rank_bundle(F::FilteredBundle) -> Int

Return the total rank of the filtered bundle.
"""
function rank_bundle(F::FilteredBundle)
  sum(rank_bundle(p) for p in F.pieces; init=0)
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
function _root_heights(mdt::MarkedDynkinType)
  DT = _ambient_type(mdt)
  R = rank(DT)
  RS = RootSystem(DT)
  Marked = marked_nodes(mdt)
  unmarked = unmarked_nodes(mdt)

  nonpar_roots = positive_nonparabolic_roots(mdt)

  # Group roots by nonparabolic height
  height_groups = Dict{Int,Vector{typeof(first(nonpar_roots))}}()
  for α in nonpar_roots
    h = _nonparabolic_height(coefficients(α), Marked)
    push!(get!(height_groups, h, typeof(α)[]), α)
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

The tangent bundle ``\\mathrm{T}_{G/P}`` with its natural filtration by
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
function filtered_tangent_bundle(X::PartialFlagVariety)
  mdt = marked_dynkin_type(X)
  height_data = _root_heights(mdt)

  FilteredBundle(
    X,
    CompletelyReducibleBundle[CompletelyReducibleBundle(X, ws) for (_, ws) in height_data],
  )
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
Shared construction of ``\\wedge^k F`` and ``\\mathrm{Sym}^k F`` for a filtered
bundle ``F``: the associated graded of the induced filtration is

```math
\\bigoplus_{|\\alpha|=k} P^{\\alpha_1}(\\mathrm{gr}_1) \\otimes \\cdots \\otimes P^{\\alpha_s}(\\mathrm{gr}_s)
```

for ``P = \\wedge`` or ``\\mathrm{Sym}``, where the term of multiexponent
``\\alpha`` sits at filtration weight ``\\sum_i i \\cdot \\alpha_i`` and terms of
equal weight form a single graded piece.
"""
function _graded_power(power, F::FilteredBundle, k::Integer)
  k = Int(k)
  k < 0 && return FilteredBundle(F.variety, CompletelyReducibleBundle[])
  k == 0 && return FilteredBundle(F.variety, [structure_sheaf(F.variety)])
  k == 1 && return F

  pieces = graded_pieces(F)
  s = length(pieces)
  weight_terms = Dict{Int,Vector{IrrepLevi}}()

  for α in multiexponents(s, k)
    factors = CompletelyReducibleBundle[]
    for i in 1:s
      f = power(pieces[i], α[i])
      iszero(f) && (empty!(factors); break)
      push!(factors, f)
    end
    isempty(factors) && continue
    term = reduce(tensor_product, factors)
    fw = sum(i * α[i] for i in 1:s)
    append!(get!(weight_terms, fw, IrrepLevi[]), components(term))
  end

  FilteredBundle(
    F.variety,
    CompletelyReducibleBundle[
      CompletelyReducibleBundle(F.variety, weight_terms[fw]) for
      fw in sort!(collect(keys(weight_terms)))
    ],
  )
end

"""
    exterior_power(F::FilteredBundle, k::Integer) -> FilteredBundle

The ``k``-th exterior power ``\\bigwedge^k F`` of a filtered bundle, equipped
with the induced filtration (see [`_graded_power`](@ref)).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> rank_bundle(exterior_power(F, 2)) == binomial(dimension(X), 2)
true
```
"""
exterior_power(F::FilteredBundle, k::Integer) = _graded_power(exterior_power, F, k)

"""
    symmetric_power(F::FilteredBundle, k::Integer) -> FilteredBundle

The ``k``-th symmetric power ``\\mathrm{Sym}^k F`` of a filtered bundle,
equipped with the induced filtration (see [`_graded_power`](@ref)).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> F = filtered_tangent_bundle(X);

julia> rank_bundle(symmetric_power(F, 2)) == binomial(dimension(X) + 1, 2)
true
```
"""
symmetric_power(F::FilteredBundle, k::Integer) = _graded_power(symmetric_power, F, k)

"""
    dual(F::FilteredBundle) -> FilteredBundle

The dual of a filtered bundle. The filtration is reversed: if ``F`` has
pieces ``\\mathrm{gr}_1, \\ldots, \\mathrm{gr}_s`` (bottom to top), then
``F^\\vee`` has pieces ``\\mathrm{gr}_s^\\vee, \\ldots, \\mathrm{gr}_1^\\vee``
(the dual reverses the filtration order).
"""
function dual(F::FilteredBundle)
  FilteredBundle(F.variety, [dual(p) for p in reverse(F.pieces)])
end

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, F::FilteredBundle)
  n = n_filtration_steps(F)
  r = rank_bundle(F)
  print(io, "FilteredBundle(rank $r, $n layer(s))")
end
