# ═══════════════════════════════════════════════════════════════════════════════
#  Projections — pullback and pushforward along q: G/P_J → G/P_I
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
    _projection_nodes(X::PartialFlagVariety, D::PartialFlagVariety) -> Tuple{Vararg{Int}}

Check that the marked nodes ``I`` of `X` are contained in the marked nodes ``J``
of `D`, so that ``\\mathrm{P}_J \\subseteq \\mathrm{P}_I`` and the identity of
``\\mathrm{G}`` induces a projection ``q \\colon D \\to X``, and return the
contracted nodes ``J \\setminus I``.
"""
function _projection_nodes(X::PartialFlagVariety, D::PartialFlagVariety)
  dynkin_type(X) === dynkin_type(D) ||
    throw(ArgumentError("$D and $X have different ambient Dynkin types."))
  I, J = marked_nodes(X), marked_nodes(D)
  issubset(I, J) || throw(
    ArgumentError(
      "$D → $X is not a projection: the marked nodes $I of the target are not contained in the marked nodes $J of the source."
    ),
  )
  Tuple(j for j in J if !(j in I))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Pullback
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _graded_branching!(pieces::Dict{Int, Vector{IrrepLevi}}, mdt_D::MarkedDynkinType, rep::IrrepLevi)

Restrict the fibre ``\\mathrm{V}_{\\mathrm{L}_I}(\\lambda)`` of the irreducible
`rep` on ``\\mathrm{G}/\\mathrm{P}_I`` to ``\\mathrm{L}_J``, and add the
resulting irreducibles to `pieces`, keyed by the grading ``g`` of
[`pullback`](@ref).

Write a weight of the fibre as ``\\lambda - \\sum_k d_k \\alpha_{o(k)}``, where
``o`` sends the ``k``-th node of the ``\\mathrm{L}_I`` diagram to its ambient
node and the ``d_k`` are non-negative.  The coefficients at the contracted nodes
are constant on each ``\\mathrm{L}_J``-irreducible, since two weights lie in the
same one only if they differ by a root of ``\\mathrm{L}_J``.  Bucketing by those
coefficients therefore separates the irreducibles enough for
`Semisimple.character_from_weights` to peel each bucket on its own, and within a
bucket the ``\\mathrm{L}_J``-weight determines the ambient weight, so the ambient
lift survives the peeling.
"""
function _graded_branching!(
  pieces::Dict{Int,Vector{IrrepLevi}}, mdt_D::MarkedDynkinType, rep::IrrepLevi
)
  mdt_X = marked_dynkin_type(rep)
  C = cartan_matrix(dynkin_type(mdt_X))
  Cinv = cartan_matrix_inverse(levi_type(mdt_X))
  λ = coefficients(p_dominant_weight(rep))
  ss = semisimple_part(rep)
  ord = levi_permutation(mdt_X)     # node of the L_I diagram ↦ ambient node
  ord_D = levi_permutation(mdt_D)   # node of the L_J diagram ↦ ambient node
  contracted = Tuple(k for k in eachindex(ord) if !(ord[k] in ord_D))

  # One vector per weight of the fibre, so the cost scales with the number of
  # weights of V_{L_I}(λ) and the bookkeeping dominates: 28 MiB for a
  # 3248-weight fibre on the Cayley plane, of which only 1 MiB is the weight
  # system itself.  Keeping the per-weight data in tuples would cut the constant
  # factor; a branching rule producing highest weights directly would remove the
  # |W_{L_J}| factor altogether.
  buckets = Dict{Vector{Int},Vector{Pair{Vector{Int},BigInt}}}()
  for (ν, mult) in freudenthal_formula(ss)
    d = Int.(Cinv * (coefficients(ss) - ν))   # the d_k, since ss - ν = Σ_k d_k α_k
    w = Int[λ[i] - sum(d[k] * C[i, ord[k]] for k in eachindex(ord)) for i in eachindex(λ)]
    push!(
      get!(buckets, Int[d[k] for k in contracted], Pair{Vector{Int},BigInt}[]), w => mult
    )
  end

  LT_D = levi_type(mdt_D)
  for (key, entries) in buckets
    reps = get!(pieces, sum(key), IrrepLevi[])
    if LT_D === nothing
      # D = G/B: the Levi is a torus, so each bucket is a single fibre weight.
      append!(reps, (IrrepLevi(mdt_D, w) for (w, m) in entries for _ in 1:m))
      continue
    end
    r = rank(LT_D)
    coords(w) = SVector{r,Int}(Int[w[j] for j in ord_D])
    χ = character_from_weights(
      LT_D, Dict{SVector{r,Int},BigInt}(coords(w) => m for (w, m) in entries)
    )
    lift = Dict(coords(w) => w for (w, _) in entries)
    append!(reps, (IrrepLevi(mdt_D, lift[coefficients(ω)]) for (ω, m) in χ for _ in 1:m))
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
[`_graded_branching!`](@ref).  Then
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
  # Nothing is contracted: q is the identity, up to the display name of D.
  isempty(_projection_nodes(variety(E), D)) &&
    return FilteredBundle(D, [CompletelyReducibleBundle(D, components(E))])

  mdt_D = marked_dynkin_type(D)
  pieces = Dict{Int,Vector{IrrepLevi}}()
  for rep in components(E)
    # A non-dominant weight is the zero bundle, by the `fiber_dimension` convention.
    iszero(fiber_dimension(rep)) && continue
    _graded_branching!(pieces, mdt_D, rep)
  end

  FilteredBundle(
    D,
    CompletelyReducibleBundle[
      CompletelyReducibleBundle(D, pieces[g]) for
      g in sort!(collect(keys(pieces)); rev=true)
    ],
  )
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Pushforward
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _fold_dominant_in(μ::WeightLatticeElem, nodes) -> Tuple{WeightLatticeElem, Int}

Fold the ambient weight `μ` into the dominant chamber of the root subsystem
spanned by the simple roots in `nodes`, and return it together with the length of
the element of ``\\mathrm{W}_S`` used.  For `nodes` all of them this is
`Semisimple.conjugate_dominant_weight_with_length`.
"""
function _fold_dominant_in(μ::WeightLatticeElem, nodes)
  len = 0
  while (k = findfirst(i -> coefficients(μ)[i] < 0, nodes)) !== nothing
    μ = reflect(μ, nodes[k])
    len += 1
  end
  (μ, len)
end

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
dominant for ``\\mathrm{L}_I``.  If ``\\lambda+\\rho`` is singular for a root of
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
  _projection_nodes(X, variety(E))
  mdt_X = marked_dynkin_type(X)
  S = unmarked_nodes(mdt_X)
  ρ = weyl_vector(dynkin_type(X))
  d = dimension(variety(E)) - dimension(X)

  pieces = [IrrepLevi[] for _ in 0:d]
  for rep in components(E)
    μ, len = _fold_dominant_in(p_dominant_weight(rep) + ρ, S)
    any(s -> iszero(coefficients(μ)[s]), S) && continue   # singular for L_I
    len <= d || continue
    push!(pieces[len + 1], IrrepLevi(mdt_X, μ - ρ))
  end

  Cohomology{CompletelyReducibleBundle}(
    CompletelyReducibleBundle[CompletelyReducibleBundle(X, p) for p in pieces], d
  )
end

function Base.show(io::IO, Rq::Cohomology{CompletelyReducibleBundle})
  parts = [
    "R$(_superscript(i)) = $(Rq[i])" for i in 0:lastindex(Rq) if
                                         !isempty(components(Rq[i]))
  ]
  print(io, isempty(parts) ? "R* = 0" : join(parts, "\n"))
end
