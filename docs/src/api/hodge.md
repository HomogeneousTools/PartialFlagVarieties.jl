# Hodge & Hochschild

Hodge numbers, twisted Hodge numbers, and Hochschild cohomology.

See the [mathematical background](../math.md#Hochschild-cohomology-of-zero-loci) for a detailed
description of the algorithm, including why both twisted Hodge matrices ``M_1`` and ``M_2``
are needed to resolve symbolic variables, how Euler characteristic constraints supplement
the Serre cross-constraints, and the distinction between ample and nef anticanonical bundles
in the Fano vanishing check.

There are two rather different regimes:

- for homogeneous varieties ``X = \mathrm{G}/\mathrm{P}``, Hodge-theoretic output is integer-valued;
- for zero loci ``Z``, the package may return symbolic affine expressions when
  the long exact sequences do not force all connecting-map ranks.

## Hodge numbers

Both `hodge_numbers(X)` and `twisted_hodge_numbers(X, L)` return an
`OffsetArrays.OffsetMatrix` with axes `0:dimension(X)` in both directions.
Thus `H[p, q]` is the Hodge number with indices `(p, q)`, without adding one.
The aliases `hodge_numbers_symbolic` and `hodge_numbers_les` use the same convention.

This changes the previous one-based API: replace `H[p+1, q+1]` by `H[p, q]`.
Use `parent(H)` to access the underlying ordinary one-based matrix.

```@docs
hodge_numbers
twisted_hodge_numbers
print_hodge_diamond
```

## Hochschild cohomology

The Hochschild–Kostant–Rosenberg decomposition gives

```math
\mathrm{HH}^n(X) = \bigoplus_{p+q=n} \mathrm{H}^q(X, \bigwedge\nolimits^p \mathrm{T}_X)
```

The result is stored in a [`PolyvectorParallelogram`](@ref), which uses
0-based access: `P[p, q] = h^q(X, \bigwedge^p T_X)` for individual entries and
`P[n] = \dim \mathrm{HH}^n(X)` for the total dimension along the anti-diagonal.

```@docs
hochschild_cohomology
PolyvectorParallelogram
```

## Internals

The conormal restriction backends and the Lefschetz refinement used by
[`hodge_numbers`](@ref).

```@docs
PartialFlagVarieties.GradedConormal
PartialFlagVarieties.FilteredConormal
PartialFlagVarieties._lefschetz_inject!
```
