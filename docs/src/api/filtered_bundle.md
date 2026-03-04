# Filtered Bundles

Equivariant bundles equipped with a filtration by subbundles.

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
