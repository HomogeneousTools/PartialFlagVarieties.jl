# Spectral Sequences

A [`FilteredBundle`](@ref) ``\mathcal{F}`` on ``X = \mathrm{G}/\mathrm{P}`` determines a
spectral sequence

```math
\mathrm{E}_1^{p,q} \Longrightarrow \mathrm{H}^{p+q}(X, \mathcal{F}),
```

where ``\mathrm{E}_1^{p,q}`` is the degree-``(p+q)`` cohomology of the
``q``-th graded piece of ``\mathcal{F}``, counted from the top of the
filtration, computed by Borel–Weil–Bott.  The page is supported on the
horizontal band ``0 \le q \le s - 1``, ``0 \le p + q \le \dim X`` (with
``s`` the number of graded pieces) — not on the first quadrant, since ``p``
may be negative — and every differential
``d_r \colon \mathrm{E}_r^{p,q} \to \mathrm{E}_r^{p+1-r,\, q+r}`` raises
the total degree by exactly one and strictly raises ``q``. On a
cominuscule ``\mathrm{G}/\mathrm{P}`` (abelian nilradical) the tangent bundle is completely
reducible and nothing can degenerate; in general the differentials may be
nonzero, and computing with the associated graded alone silently assumes
degeneration. This module makes that step honest: all differentials are
``\mathrm{G}``-equivariant, so the spectral sequence splits into isotypical components,
and both the degeneracy test and the cohomology of the abutment
(`cohomology(F::FilteredBundle)`, with symbolic entries where differentials
remain possible) work per component.

The same mechanism drives the Hodge number computations for zero loci in
non-cominuscule ambient varieties.

## Types

```@docs
SpectralSequence
```

## Functions

```@docs
spectral_sequence
E1_page
isotypical_components
does_E1_degenerate
```
