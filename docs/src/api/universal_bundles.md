# Universal Bundles

Tautological (universal) bundles and spinor bundles on classical varieties.

The geometric meaning is straightforward on generalized Grassmannians and
quadrics. On multi-step type-A flag varieties, however,
[`tautological_bundles`](@ref) and [`universal_subbundles`](@ref) should be read as
**convenient equivariant building blocks** attached to the marked nodes, not as
literal filtered subbundles or quotients inside a full tautological flag.

## Grassmannians

```@docs
universal_subbundle
universal_quotient_bundle
S
Q
```

## Spinor bundles

```@docs
spinor_bundle
```

## Partial flag varieties

!!! warning
    For multi-step flags these functions return completely reducible bundles
    attached to the marked nodes. They agree with the usual tautological bundle
    only in the one-step Grassmannian case.

```@docs
tautological_bundles
universal_subbundles
```
