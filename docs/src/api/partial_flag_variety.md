# PartialFlagVariety

The primary user-facing type for partial flag varieties ``G/P``.

A `PartialFlagVariety` wraps a [`MarkedDynkinType`](@ref) and provides the
high-level API for computing invariants of ``G/P``: dimension, Picard rank,
Euler characteristic, Betti numbers, and classification predicates.

For the package conventions behind the constructors and invariants on this
page, see [Conventions & Notation](../conventions.md).

## Constructors

If a variety already has a standard name, prefer the
[named constructors](constructions.md) such as `Gr`, `quadric`, or
`cayley_plane`.

For an arbitrary partial flag variety, use `partial_flag_variety` directly with
either a Dynkin type or a Dynkin-type string:

```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{3}, (1, 3));

julia> picard_rank(X)
2

julia> dimension(partial_flag_variety("D5", 5))
10

julia> marked_nodes(PartialFlagVariety("A3", [2]))
(2,)
```

The two-argument form `PartialFlagVariety("A3", [2])` is just a convenience
alias for `partial_flag_variety("A3", [2])`.

!!! note
    The one-argument constructor `PartialFlagVariety("31")` is different: it
    decodes a ZeroLocus62 label rather than parsing a Dynkin-type string. See
    [Labels](labels.md).

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
