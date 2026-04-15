# Hodge & Hochschild

Hodge numbers, twisted Hodge numbers, and Hochschild cohomology.

See the [mathematical background](../math.md#Hochschild-cohomology-of-zero-loci) for a detailed
description of the algorithm, including why both twisted Hodge matrices ``M_1`` and ``M_2``
are needed to resolve symbolic variables, how Euler characteristic constraints supplement
the Serre cross-constraints, and the distinction between ample and nef anticanonical bundles
in the Fano vanishing check.

## Hodge numbers

```@docs
hodge_numbers
twisted_hodge_numbers
print_hodge_diamond
```

## Hochschild cohomology

The Hochschild–Kostant–Rosenberg decomposition gives

```math
\mathrm{HH}^n(X) = \bigoplus_{p+q=n} H^q(X, \bigwedge^p T_X)
```

```@docs
hochschild_cohomology
hochschild_dimension
PolyvectorParallelogram
```
