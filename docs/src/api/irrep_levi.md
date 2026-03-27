# IrrepLevi

Irreducible representations of the Levi subgroup ``L`` of the parabolic
``P \subseteq G``.

## Mathematical background

The Levi factor of a parabolic subgroup decomposes as
``L = \operatorname{Z}(L)^\circ \cdot [L,L]``, where ``\operatorname{Z}(L)^\circ`` is the connected center
(a torus of rank equal to the number of marked nodes ``|I|``) and ``[L,L]``
is the semisimple part (whose Dynkin diagram is the sub-diagram on the
unmarked nodes).

Every irreducible ``L``-representation therefore splits as

```math
V = \chi \otimes W
```

where ``\chi`` is a **character** of the center ``\operatorname{Z}(L)^\circ`` and ``W`` is
an **irreducible** representation of ``[L,L]`` with highest weight ``\mu``.

An `IrrepLevi` stores this decomposition: the central part is a vector of
integers (the character coordinates times the
[`central_scaling_factor`](@ref)), and the semisimple part is a
`WeightLatticeElem` for ``[L,L]``.

## Design notes

!!! note "Why `IrrepLevi` is an internal type"
    Users construct equivariant bundles via the public API:
    [`line_bundle`](@ref), [`tangent_bundle`](@ref), [`exterior_power`](@ref),
    etc. These functions build [`CompletelyReducibleBundle`](@ref) objects
    whose internal decomposition into `IrrepLevi` components happens
    automatically.

    The situations where a user encounters `IrrepLevi` directly are:
    - Inspecting [`components`](@ref)`(E)` to see the irreducible summands of
      a bundle
    - Advanced computations that need the central/semisimple coordinates

!!! note "Caching of tensor products"
    Tensor products of `IrrepLevi` decompose via Lie.jl character arithmetic,
    which involves computing tensor product branching rules from Weyl character
    data. This is the hottest path in the package and is **LRU-cached**: see
    [`PartialFlagVarieties.configure_caches!`](@ref) for tuning the cache budget.

## Type

```@docs
IrrepLevi
```

## Accessors

```@docs
central_part
semisimple_part
p_dominant_weight
fiber_dimension
to_ambient_weight
```

## Monoidal operations

```@docs
tensor_product(::IrrepLevi, ::IrrepLevi)
dual(::IrrepLevi)
exterior_power(::IrrepLevi, ::Integer)
symmetric_power(::IrrepLevi, ::Integer)
```
