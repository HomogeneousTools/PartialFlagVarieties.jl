# Filtered Bundles

Equivariant bundles equipped with a filtration by subbundles.

The primary source of filtered bundles is the **tangent bundle** of ``\mathrm{G}/\mathrm{P}``.
As a ``\mathrm{P}``-representation, the tangent space ``\mathfrak{g}/\mathfrak{p}``
has a natural **root height filtration**: the roots ``\alpha \in \Phi^+ \setminus
\Phi^+_{\mathrm{L}}`` are stratified by height (the sum of their simple root coefficients),
and this stratification yields a ``\mathrm{P}``-equivariant filtration

```math
0 = \mathcal{F}_0 \subset \mathcal{F}_1 \subset \cdots \subset \mathcal{F}_s = \mathrm{T}_{\mathrm{G}/\mathrm{P}}
```

whose graded pieces ``\mathrm{gr}_k = \mathcal{F}_k / \mathcal{F}_{k-1}`` are completely reducible
``\mathrm{L}``-representations. This filtration is key to computing
``\mathrm{H}^q(X, \Omega^p_X)`` via the associated spectral sequence, and to the
conormal filtration used in zero-locus Hodge number computations.

## API

```@docs
FilteredBundle
rank(::FilteredBundle)
graded_pieces
n_filtration_steps
total_bundle
filtered_tangent_bundle
filtered_cotangent_bundle
```

## Operations

Filtered bundles support exterior powers, symmetric powers, determinants,
duals, twists, and tensor products with either completely reducible or filtered
bundles. These operations preserve and refine the filtration structure.
[`external_tensor_product`](@ref) also accepts filtered bundles on different
partial flag varieties and returns the induced filtration on their product.

For a filtered bundle ``\mathcal{F}`` with graded pieces ``\mathrm{gr}_1, \ldots, \mathrm{gr}_s``,
the ``k``-th exterior power ``\bigwedge\nolimits^k \mathcal{F}`` has an induced filtration whose
associated graded is

```math
\bigoplus_{|\alpha|=k} \bigwedge\nolimits^{\alpha_1} \mathrm{gr}_1 \otimes \cdots \otimes \bigwedge\nolimits^{\alpha_s} \mathrm{gr}_s
```

ordered by the filtration weight ``\sum_i i \cdot \alpha_i``.

```@docs
exterior_power(::FilteredBundle, ::Integer)
symmetric_power(::FilteredBundle, ::Integer)
dual(::FilteredBundle)
tensor_product(::FilteredBundle, ::CompletelyReducibleBundle)
tensor_product(::FilteredBundle, ::FilteredBundle)
det(::FilteredBundle)
twist(::FilteredBundle, ::Integer, ::Integer)
Base.iszero(::FilteredBundle)
```

## Internals

```@docs
PartialFlagVarieties._graded_power
```
