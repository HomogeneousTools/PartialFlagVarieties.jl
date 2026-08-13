# Projections

Pullback and pushforward along a projection between partial flag varieties.

Marked nodes ``I \subseteq J`` give ``\mathrm{P}_J \subseteq \mathrm{P}_I``, hence a
projection

```math
q \colon D = \mathrm{G}/\mathrm{P}_J \longrightarrow \mathrm{G}/\mathrm{P}_I = X
```

whose fibre ``\mathrm{P}_I/\mathrm{P}_J \cong \mathrm{L}_I/(\mathrm{L}_I \cap \mathrm{P}_J)``
is again a partial flag variety, for the Levi factor ``\mathrm{L}_I`` in place of
``\mathrm{G}``. Both directions along ``q`` stay inside the world of equivariant
bundles, and each is a familiar construction relative to ``\mathrm{L}_I``.

## Pullback

The unipotent radical of ``\mathrm{P}_J`` acts nontrivially on the fibre of
``q^*\mathcal{E}``, so ``q^*\mathcal{E}`` is in general only a
[`FilteredBundle`](@ref). Its associated graded is completely reducible, being
the restriction of the fibre ``\mathrm{V}_{\mathrm{L}_I}(\lambda)`` to the
reductive subgroup ``\mathrm{L}_J``. This is Levi branching, where the Cartan
subgroup is shared, so it is pure weight bookkeeping: list the weights of the
fibre, lift them to the ambient weight lattice, and collect them into
``\mathrm{L}_J``-irreducibles.

The graded pieces are ordered from subbundle to quotient, which is the order
[`graded_pieces`](@ref) expects. On ``\mathrm{Fl}(k, l; n)`` the pullback of the
tautological subbundle of ``\mathrm{Gr}(l, n)`` is the tautological filtration

```math
0 \to \mathcal{U}_k \to \mathcal{U}_l \to \mathcal{U}_l/\mathcal{U}_k \to 0,
```

with ``\mathcal{U}_k`` first.

```@docs
pullback
```

## Pushforward

``\mathrm{R}q_*`` is the Borel–Weil–Bott theorem for ``\mathrm{L}_I`` rather than
for ``\mathrm{G}``: the same ``\rho``-shifted fold into a dominant chamber, but
reflecting only in the nodes unmarked in ``I``. As in the absolute case, an
irreducible pushes forward to a single shifted irreducible or to zero, so no
branching is involved. Taking ``I = \emptyset``, where ``X`` is a point, recovers
[`cohomology`](@ref).

Since ``\mathrm{R}q_*`` of an irreducible sits in a single degree ``d``, the Leray
spectral sequence degenerates and

```math
\mathrm{H}^i(D, \mathcal{E}) = \mathrm{H}^{i-d}(X, \mathrm{R}^dq_*\mathcal{E}).
```

```@docs
pushforward
```

## Internals

```@docs
PartialFlagVarieties._check_projection
PartialFlagVarieties._graded_branching
PartialFlagVarieties._fibre_buckets
```
