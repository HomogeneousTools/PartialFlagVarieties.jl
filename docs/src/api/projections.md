# Projections

Pullback along a projection between partial flag varieties.

Marked nodes ``I \subseteq J`` give ``\mathrm{P}_J \subseteq \mathrm{P}_I``, hence a
projection

```math
q \colon D = \mathrm{G}/\mathrm{P}_J \longrightarrow \mathrm{G}/\mathrm{P}_I = X
```

whose fibre ``\mathrm{P}_I/\mathrm{P}_J \cong \mathrm{L}_I/(\mathrm{L}_I \cap \mathrm{P}_J)``
is again a partial flag variety, for the Levi factor ``\mathrm{L}_I`` in place of
``\mathrm{G}``.

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

## Internals

```@docs
PartialFlagVarieties._projection_nodes
PartialFlagVarieties._graded_branching!
```
