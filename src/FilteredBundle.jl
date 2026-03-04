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
function tensor_product(F::FilteredBundle{MDT}, E::CompletelyReducibleBundle{MDT}) where {MDT}
  new_pieces = [tensor_product(p, E) for p in F.pieces]
  FilteredBundle(F.variety, new_pieces)
end

function tensor_product(E::CompletelyReducibleBundle{MDT}, F::FilteredBundle{MDT}) where {MDT}
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

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, F::FilteredBundle{MDT}) where {MDT}
  n = n_filtration_steps(F)
  r = rank_bundle(F)
  print(io, "FilteredBundle(rank $r, $n layer(s))")
end
