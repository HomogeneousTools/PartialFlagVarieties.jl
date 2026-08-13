# ═══════════════════════════════════════════════════════════════════════════════
#  Projections: pullback and pushforward along q: G/P_J → G/P_I
#
#  For I ⊆ J we have P_J ⊆ P_I, hence a projection
#
#      q: D = G/P_J ⟶ G/P_I = X
#
#  whose fibre P_I/P_J ≅ L_I/(L_I ∩ P_J) is itself a partial flag variety, for
#  the Levi factor L_I instead of for G.  Both directions along q stay inside
#  the world of equivariant bundles:
#
#   * q^*E is in general only filtered, because the unipotent radical of P_J
#     acts nontrivially on the fibre.  The associated graded is the branching of
#     the fibre from L_I to L_J, which is completely reducible since L_J is
#     reductive.
#   * Rq_* is the Borel–Weil–Bott theorem for L_I rather than for G: the same
#     ρ-shifted fold into a dominant chamber, but reflecting only in the nodes
#     unmarked in I.  Taking I = ∅ recovers `cohomology`.
# ═══════════════════════════════════════════════════════════════════════════════

export pullback, pushforward

"""
    _check_projection(X::PartialFlagVariety, D::PartialFlagVariety)

Check that the marked nodes ``I`` of `X` are contained in the marked nodes ``J``
of `D`, so that ``\\mathrm{P}_J \\subseteq \\mathrm{P}_I`` and the identity of
``\\mathrm{G}`` induces a projection ``q \\colon D \\to X``.
"""
function _check_projection(X::PartialFlagVariety, D::PartialFlagVariety)
  dynkin_type(X) === dynkin_type(D) ||
    throw(ArgumentError("$D and $X have different ambient Dynkin types."))
  I, J = marked_nodes(X), marked_nodes(D)
  issubset(I, J) || throw(
    ArgumentError(
      "$D → $X is not a projection: the marked nodes $I of the target are not contained in the marked nodes $J of the source."
    ),
  )
  return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Pullback
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _graded_branching(mdt_D::MarkedDynkinType, rep::IrrepLevi) -> Dict{Int, Vector{IrrepLevi}}

Restrict the fibre ``\\mathrm{V}_{\\mathrm{L}_I}(\\lambda)`` of the irreducible
`rep` on ``\\mathrm{G}/\\mathrm{P}_I`` to ``\\mathrm{L}_J``, and return the
resulting irreducibles keyed by the grading ``g`` of [`pullback`](@ref).

Bucket the weights of the fibre, then peel each bucket into irreducibles: see
[`_fibre_buckets`](@ref) for why bucketing first is what makes the peeling work.
"""
function _graded_branching(mdt_D::MarkedDynkinType, rep::IrrepLevi)
  mdt_X = marked_dynkin_type(rep)
  buckets = _fibre_buckets(
    dynkin_type(mdt_X), levi_type(mdt_X), rep, unmarked_nodes(mdt_D)
  )

  pieces = Dict{Int,Vector{IrrepLevi}}()
  for (key, entries) in buckets
    reps = get!(pieces, sum(key), IrrepLevi[])
    _peel!(reps, levi_type(mdt_D), mdt_D, entries)
  end
  pieces
end

"""
    _fibre_buckets(::Type{DT}, ::Type{LT_I}, rep::IrrepLevi, kept) -> Dict

Group the weights of the fibre of `rep` into buckets, each holding the ambient
weights it contains with their multiplicities. `kept` is the set of ambient
nodes surviving in ``\\mathrm{L}_J``, i.e. the nodes unmarked in ``J``.

Every weight of the fibre is
``\\lambda - \\sum_k d_k \\alpha_{o(k)}`` for non-negative integers ``d_k``,
where ``o`` sends the ``k``-th node of the ``\\mathrm{L}_I`` diagram to its
ambient node; the bucket key is the vector of those ``d_k`` at the *contracted*
nodes, the ones in ``J \\setminus I``.

That key is constant on each ``\\mathrm{L}_J``-irreducible, because two weights
lie in the same one only if they differ by a root of ``\\mathrm{L}_J``. So the
buckets separate the irreducibles enough for `Semisimple.character_from_weights`
to peel each on its own, and within a bucket the ``\\mathrm{L}_J``-weight
determines the ambient weight, so the ambient lift survives the peeling.

The Dynkin types are taken as type parameters rather than read off the marked
Dynkin type, where they are runtime values: otherwise the Cartan matrices infer
as `Any` and this loop runs through dynamic dispatch.

The whole weight system is materialised, one vector per weight, so the cost
scales with the number of weights of ``\\mathrm{V}_{\\mathrm{L}_I}(\\lambda)``
and the bookkeeping dominates. A branching rule producing highest weights
directly would remove the ``|\\mathrm{W}_{\\mathrm{L}_J}|`` factor altogether.
"""
function _fibre_buckets(
  ::Type{DT}, ::Type{LT_I}, rep::IrrepLevi, kept
) where {DT<:DynkinType,LT_I<:DynkinType}
  C = cartan_matrix(DT)
  Cinv = cartan_matrix_inverse(LT_I)
  λ = coefficients(p_dominant_weight(rep))
  ss = semisimple_part(rep)
  ord = levi_permutation(marked_dynkin_type(rep))
  contracted = Tuple(k for k in eachindex(ord) if !(ord[k] in kept))

  buckets = Dict{Vector{Int},Vector{Pair{Vector{Int},BigInt}}}()
  for (ν, mult) in freudenthal_formula(ss)
    d = Int.(Cinv * (coefficients(ss) - ν))
    w = Int[λ[i] - sum(d[k] * C[i, ord[k]] for k in eachindex(ord)) for i in eachindex(λ)]
    key = Int[d[k] for k in contracted]
    push!(get!(() -> Pair{Vector{Int},BigInt}[], buckets, key), w => mult)
  end
  buckets
end

# D = G/B: the Levi is a torus, so each bucket is a single fibre weight.
function _peel!(reps, ::Nothing, mdt_D::MarkedDynkinType, entries)
  append!(reps, (IrrepLevi(mdt_D, w) for (w, m) in entries for _ in 1:Int(m)))
end

# Otherwise peel the bucket into L_J-irreducibles, looking each ambient lift
# back up by the L_J-weight that identifies it.
function _peel!(
  reps, ::Type{LT_J}, mdt_D::MarkedDynkinType, entries
) where {LT_J<:DynkinType}
  ord_D = levi_permutation(mdt_D)
  r = rank(LT_J)
  mults = Dict{SVector{r,Int},BigInt}()
  lift = Dict{SVector{r,Int},Vector{Int}}()
  for (w, m) in entries
    c = SVector{r,Int}(ntuple(j -> w[ord_D[j]], r))
    mults[c] = get(mults, c, zero(BigInt)) + m   # injective here; accumulate anyway
    lift[c] = w
  end
  χ = character_from_weights(LT_J, mults)
  for (ω, m) in χ
    # `character_from_weights` returns a virtual character in general; a negative
    # multiplicity here would mean the fibre was not a genuine representation,
    # and `1:(-2)` would drop it without trace.
    m > 0 || throw(ArgumentError("Branching produced multiplicity $m for $ω."))
    append!(reps, Iterators.repeated(IrrepLevi(mdt_D, lift[coefficients(ω)]), Int(m)))
  end
end

"""
    pullback(D::PartialFlagVariety, E::CompletelyReducibleBundle) -> FilteredBundle

Pull ``\\mathcal{E}`` back along the projection
``q \\colon D = \\mathrm{G}/\\mathrm{P}_J \\to \\mathrm{G}/\\mathrm{P}_I = X``,
where `E` lives on ``X`` and the marked nodes of `X` are contained in those of
`D`.

The result is only a [`FilteredBundle`](@ref): the unipotent radical of
``\\mathrm{P}_J`` acts nontrivially on the fibre of ``q^*\\mathcal{E}``, so
``q^*\\mathcal{E}`` does not split, but its associated graded does, being the
restriction of the fibre from ``\\mathrm{L}_I`` to the reductive subgroup
``\\mathrm{L}_J``.

# Filtration order

Write a weight of the fibre as ``\\lambda - \\sum_k d_k \\alpha_{o(k)}``, as in
[`_graded_branching`](@ref).  Then
``g = \\sum_{o(k) \\in J \\setminus I} d_k`` is constant on each graded piece, and
the pieces are returned by *descending* ``g``.  This is the order
[`graded_pieces`](@ref) expects: the piece of maximal ``g`` is the subbundle, and
the piece containing the highest weight ``\\lambda``, where ``g = 0``, is the top
quotient.  The natural guess is the other way around.

# Examples

On ``\\mathrm{Fl}(1, 3; 6)`` the pullback of the tautological subbundle of
``\\mathrm{Gr}(3, 6)`` is the tautological filtration
``0 \\to \\mathcal{U}_1 \\to \\mathcal{U}_3 \\to \\mathcal{U}_3/\\mathcal{U}_1 \\to 0``:

```jldoctest
julia> using PartialFlagVarieties

julia> D = flag_variety(6, [1, 3]);

julia> F = pullback(D, universal_subbundle(Gr(3, 6)));

julia> graded_pieces(F) == tautological_bundles(D)
true

julia> rank_bundle.(graded_pieces(F))
2-element Vector{Int64}:
 1
 2
```

Pulling back a line bundle keeps it a line bundle, of degree zero at the
contracted nodes:

```jldoctest
julia> using PartialFlagVarieties

julia> F = pullback(flag_variety(4, [1, 2]), line_bundle(Gr(2, 4), 3));

julia> picard_degrees(total_bundle(F))
2-element Vector{Int64}:
 0
 3
```
"""
function pullback(D::PartialFlagVariety, E::CompletelyReducibleBundle)
  _check_projection(variety(E), D)
  mdt_D = marked_dynkin_type(D)

  # I == J: q is the identity, and the only case where the source Levi can be a
  # torus, which the branching below has no diagram to work with.  Drop the
  # fibreless summands as the general path does, so the two agree.
  if marked_dynkin_type(variety(E)) == mdt_D
    reps = IrrepLevi[rep for rep in components(E) if has_fiber(rep)]
    return FilteredBundle(
      D,
      isempty(reps) ? CompletelyReducibleBundle[] : [CompletelyReducibleBundle(D, reps)],
    )
  end

  pieces = Dict{Int,Vector{IrrepLevi}}()
  # Plethysms and tensor products repeat a summand once per multiplicity, so
  # branch each distinct component once and reuse it.  Walking `components(E)`
  # in order keeps the order within each graded piece.
  branched = Dict{IrrepLevi,Dict{Int,Vector{IrrepLevi}}}()
  for rep in components(E)
    # A non-dominant weight is the zero bundle, by the `fiber_dimension` convention.
    has_fiber(rep) || continue
    for (g, reps) in get!(() -> _graded_branching(mdt_D, rep), branched, rep)
      append!(get!(pieces, g, IrrepLevi[]), reps)
    end
  end

  FilteredBundle(
    D,
    CompletelyReducibleBundle[
      CompletelyReducibleBundle(D, pieces[g]) for
      g in sort!(collect(keys(pieces)); rev=true)
    ],
  )
end

"""
    pullback(D::PartialFlagVariety, F::FilteredBundle) -> FilteredBundle

Pull a filtered bundle back by pulling back each graded piece and concatenating
the results, so that the filtration of ``q^*\\mathcal{F}`` refines the pullback
of the filtration of ``\\mathcal{F}``.

Pullback is exact, so this is again a filtration of ``q^*\\mathcal{F}``, and the
sub-to-quotient order of [`graded_pieces`](@ref) is preserved: the pieces of
each `pullback(D, p)` are already ordered that way, and the pieces of `F` are
traversed in that order.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> Q = universal_quotient_bundle(OGr(2, 7));   # already filtered

julia> F = pullback(partial_flag_variety(TypeB{3}, (1, 2)), Q);

julia> rank_bundle(F) == rank_bundle(Q)
true
```

Pulling back in two steps refines the filtration of doing it in one:

```jldoctest
julia> using PartialFlagVarieties

julia> U = universal_subbundle(Gr(2, 4));

julia> F = pullback(full_flag_variety(TypeA{3}), pullback(flag_variety(4, [1, 2]), U));

julia> rank_bundle.(graded_pieces(F))
2-element Vector{Int64}:
 1
 1
```
"""
function pullback(D::PartialFlagVariety, F::FilteredBundle)
  FilteredBundle(
    D,
    reduce(
      vcat,
      (graded_pieces(pullback(D, p)) for p in graded_pieces(F));
      init=CompletelyReducibleBundle[],
    ),
  )
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Pushforward
# ═══════════════════════════════════════════════════════════════════════════════

"""
    pushforward(X::PartialFlagVariety, E::CompletelyReducibleBundle) -> Cohomology{CompletelyReducibleBundle}

Compute ``\\mathrm{R}q_*\\mathcal{E}`` along the projection
``q \\colon D = \\mathrm{G}/\\mathrm{P}_J \\to \\mathrm{G}/\\mathrm{P}_I = X``,
where `E` lives on ``D`` and the marked nodes of `X` are contained in those of
`D`.

The result is indexed by the degree of the higher direct image, so `Rq[i]` is
``\\mathrm{R}^iq_*\\mathcal{E}`` as a bundle on `X`, for
``0 \\le i \\le \\dim D - \\dim X``.

# Algorithm

The fibre of ``q`` is
``\\mathrm{P}_I/\\mathrm{P}_J \\cong \\mathrm{L}_I/(\\mathrm{L}_I \\cap \\mathrm{P}_J)``,
so this is Borel–Weil–Bott applied fibrewise for ``\\mathrm{L}_I``.  For an
irreducible summand of ambient weight ``\\lambda``, let ``S`` be the nodes
unmarked in ``I`` and find ``w \\in \\mathrm{W}_S`` making ``w(\\lambda+\\rho)``
dominant for ``\\mathrm{L}_I``, which is what
[`borel_weil_bott(λ, nodes)`](@ref) does.  If ``\\lambda+\\rho`` is singular for a root of
``S`` then ``\\mathrm{R}q_*`` of that summand vanishes; otherwise the summand
contributes the bundle of ambient weight ``w(\\lambda+\\rho) - \\rho`` in degree
``\\ell(w)``, and nothing in any other degree.

Two simplifications are built in: ``\\rho = \\rho_{\\mathrm{G}}`` may be used
instead of ``\\rho_S``, since the difference pairs to zero with every coroot in
``S`` and is therefore ``\\mathrm{W}_S``-invariant; and the resulting weight is
automatically ``\\mathrm{P}_I``-dominant.

# Examples

Contracting ``\\mathrm{Fl}(1, 2; 4) \\to \\mathrm{Gr}(2, 4)``, the structure sheaf
pushes forward to the structure sheaf, in degree zero:

```jldoctest
julia> using PartialFlagVarieties

julia> Rq = pushforward(Gr(2, 4), structure_sheaf(flag_variety(4, [1, 2])));

julia> Rq[0] == structure_sheaf(Gr(2, 4))
true

julia> Rq[1]
0
```

Taking ``I = \\emptyset``, so that ``X`` is a point, recovers
[`cohomology`](@ref): the rank of ``\\mathrm{R}^iq_*`` is ``\\dim \\mathrm{H}^i``.

```jldoctest
julia> using PartialFlagVarieties

julia> E = dual(universal_subbundle(Gr(2, 4)));

julia> Rq = pushforward(partial_flag_variety(TypeA{3}, ()), E);

julia> [rank_bundle(Rq[i]) for i in 0:4] == [Int(cohomology(E)[i]) for i in 0:4]
true
```
"""
function pushforward(X::PartialFlagVariety, E::CompletelyReducibleBundle)
  _check_projection(X, variety(E))
  mdt_X = marked_dynkin_type(X)
  S = unmarked_nodes(mdt_X)
  d = dimension(variety(E)) - dimension(X)

  pieces = [IrrepLevi[] for _ in 0:d]
  for rep in components(E)
    # A non-dominant weight is the zero bundle, by the `fiber_dimension`
    # convention; the relative singularity test below only inspects the nodes in
    # S, so it cannot see a wall at a node marked in I.
    has_fiber(rep) || continue
    result = borel_weil_bott(p_dominant_weight(rep), S)
    result === nothing && continue          # λ + ρ singular for L_I
    len, μ = result
    # `len` is in range because λ is P_J-dominant, which the `has_fiber` guard
    # above is what enforces: for such a λ the fold is a minimal-length coset
    # representative, so ℓ(w) <= |Φ_S⁺| - |Φ_T⁺| = d.
    push!(pieces[len + 1], IrrepLevi(mdt_X, μ))
  end

  Cohomology{CompletelyReducibleBundle}(
    CompletelyReducibleBundle[CompletelyReducibleBundle(X, p) for p in pieces], d
  )
end

"""
    dimensions(Rq::Cohomology{CompletelyReducibleBundle}) -> Cohomology{BigInt}

Replace each higher direct image by its rank, as [`dimensions`](@ref) does for
character-valued cohomology.

The result is indexed by the degree of the direct image and runs to the relative
dimension, so it is not cohomology of `X` despite printing with an `H`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> Rq = pushforward(Gr(2, 4), structure_sheaf(flag_variety(4, [1, 2])));

julia> dimensions(Rq)[0]
1
```
"""
function dimensions(Rq::Cohomology{CompletelyReducibleBundle})
  Cohomology{BigInt}(BigInt[rank_bundle(E) for E in Rq.entries], Rq.max_degree)
end

"""
    euler_characteristic(Rq::Cohomology{CompletelyReducibleBundle}) -> BigInt

The rank of ``\\sum_i (-1)^i [\\mathrm{R}^iq_*\\mathcal{E}]`` in the Grothendieck
group of `X`. For ``\\mathcal{E} = q^*\\mathcal{F}`` this is the rank of
``\\mathcal{F}``, since ``\\mathrm{R}q_*q^*\\mathcal{F} = \\mathcal{F}``.

This is *not* ``\\chi(X, -)`` of anything: the entries are bundles on `X`, not
cohomology groups. See the warning on [`chi`](@ref).
"""
euler_characteristic(Rq::Cohomology{CompletelyReducibleBundle}) =
  euler_characteristic(dimensions(Rq))

"""
    iszero(Rq::Cohomology{CompletelyReducibleBundle}) -> Bool

Whether every higher direct image vanishes.
"""
Base.iszero(Rq::Cohomology{CompletelyReducibleBundle}) = all(iszero, Rq.entries)

function Base.show(io::IO, Rq::Cohomology{CompletelyReducibleBundle})
  parts = [
    "R$(_superscript(i)) = $(Rq[i])" for i in 0:lastindex(Rq) if
                                         !iszero(Rq[i])
  ]
  print(io, isempty(parts) ? "R* = 0" : join(parts, "\n"))
end
