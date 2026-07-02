# Spectral Sequences

A [`FilteredBundle`](@ref) ``F`` on ``X = G/P`` determines a spectral sequence

```math
E_1 = H^*(X, \mathrm{gr}\, F) \Longrightarrow H^*(X, F)
```

whose first page is computed by Borel–Weil–Bott on the graded pieces. On a
cominuscule ``G/P`` (abelian nilradical) the tangent bundle is completely
reducible and nothing can degenerate; in general the differentials may be
nonzero, and computing with the associated graded alone silently assumes
degeneration. This module makes that step honest: all differentials are
``G``-equivariant, so the spectral sequence splits into isotypical components,
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
