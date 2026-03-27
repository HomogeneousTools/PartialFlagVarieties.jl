# PartialFlagVariety

The primary user-facing type for partial flag varieties ``G/P``.

A `PartialFlagVariety` wraps a [`MarkedDynkinType`](@ref) and provides the
high-level API for computing invariants of ``G/P``: dimension, Picard rank,
Euler characteristic, Betti numbers, and classification predicates.

## Constructors

```@docs
PartialFlagVariety
partial_flag_variety
full_flag_variety
```

## Accessing the marked data

```@docs
marked_type
```

## Topological invariants

```@docs
dimension
picard_rank
euler_characteristic(::PartialFlagVariety)
betti_numbers
anticanonical_degrees
fano_index
```

## Classification predicates

```@docs
is_generalized_grassmannian
is_full_flag_variety
is_cominuscule
is_minuscule
is_adjoint
is_coadjoint
```
