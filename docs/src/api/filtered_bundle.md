# Filtered Bundles

Equivariant bundles equipped with a filtration by subbundles.

The primary source of filtered bundles is the **tangent bundle** of ``G/P``.
As a ``P``-representation, the tangent space ``\mathfrak{g}/\mathfrak{p}``
has a natural **root height filtration**: the roots ``\alpha \in \Phi^+ \setminus
\Phi^+_L`` are stratified by height (the sum of their simple root coefficients),
and this stratification yields a ``P``-equivariant filtration

```math
0 = F_0 \subset F_1 \subset \cdots \subset F_s = \mathrm{T}_{G/P}
```

whose graded pieces ``\mathrm{gr}_k = F_k / F_{k-1}`` are completely reducible
``L``-representations. This filtration is key to computing
``H^q(X, \Omega^p_X)`` via the associated spectral sequence, and to the
conormal filtration used in zero-locus Hodge number computations.

## API

```@docs
FilteredBundle
graded_pieces
n_filtration_steps
total_bundle
filtered_tangent_bundle
```

## Operations

Filtered bundles support exterior powers, symmetric powers, duals, and
tensor products with completely reducible bundles. These operations
preserve and refine the filtration structure.

For a filtered bundle ``F`` with graded pieces ``\mathrm{gr}_1, \ldots, \mathrm{gr}_s``,
the ``k``-th exterior power ``\bigwedge^k F`` has an induced filtration whose
associated graded is

```math
\bigoplus_{|\alpha|=k} \bigwedge^{\alpha_1} \mathrm{gr}_1 \otimes \cdots \otimes \bigwedge^{\alpha_s} \mathrm{gr}_s
```

ordered by the filtration weight ``\sum_i i \cdot \alpha_i``.

```@docs
exterior_power(::FilteredBundle, ::Integer)
symmetric_power(::FilteredBundle, ::Integer)
dual(::FilteredBundle)
tensor_product(::FilteredBundle, ::CompletelyReducibleBundle)
```

## Internals

```@docs
PartialFlagVarieties._graded_power
```
