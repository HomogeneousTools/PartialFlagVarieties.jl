# CompletelyReducibleBundle

Completely reducible equivariant vector bundles on ``G/P``.

## Type

```@docs
CompletelyReducibleBundle
```

## Accessors

```@docs
variety
components
n_components
rank_bundle
```

## Standard bundles

```@docs
structure_sheaf
line_bundle
tangent_bundle
cotangent_bundle
canonical_bundle
anticanonical_bundle
det_bundle
```

## Operations

```@docs
dual(::CompletelyReducibleBundle)
tensor_product(::CompletelyReducibleBundle, ::CompletelyReducibleBundle)
direct_sum
exterior_power(::CompletelyReducibleBundle, ::Int)
symmetric_power(::CompletelyReducibleBundle, ::Int)
twist
```
