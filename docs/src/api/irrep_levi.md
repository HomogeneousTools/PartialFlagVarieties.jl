# IrrepLevi

`IrrepLevi` is the package's internal building block for irreducible Levi
representations. Most end users should **not** start here: the normal entry
points are the public bundle constructors such as [`line_bundle`](@ref),
[`tangent_bundle`](@ref), [`universal_subbundle`](@ref), and
[`CompletelyReducibleBundle`](@ref).

You will usually see `IrrepLevi` only when inspecting the summands of a bundle
via [`components`](@ref), or when doing advanced representation-theoretic work.

## Mathematical background

The Levi factor of a parabolic subgroup decomposes as
``\mathrm{L} = \operatorname{Z}(\mathrm{L})^\circ \cdot [\mathrm{L},\mathrm{L}]``, where ``\operatorname{Z}(\mathrm{L})^\circ`` is the connected center
(a torus of rank equal to the number of marked nodes ``|I|``) and ``[\mathrm{L},\mathrm{L}]``
is the semisimple part (whose Dynkin diagram is the sub-diagram on the
unmarked nodes).

Every irreducible ``\mathrm{L}``-representation therefore splits as

```math
V = \chi \otimes W
```

where ``\chi`` is a **character** of the center ``\operatorname{Z}(\mathrm{L})^\circ`` and ``W`` is
an **irreducible** representation of ``[\mathrm{L},\mathrm{L}]`` with highest weight ``\mu``.

An `IrrepLevi` stores this decomposition: the central part is a vector of
integers (the character coordinates times the
[`central_scaling_factor`](@ref)), and the semisimple part is a
`WeightLatticeElem` for ``[\mathrm{L},\mathrm{L}]``.

## Design notes

!!! note "Why `IrrepLevi` is an internal type"
    Users normally construct equivariant bundles via the public API:
    [`line_bundle`](@ref), [`tangent_bundle`](@ref), [`exterior_power`](@ref),
    [`universal_subbundle`](@ref), and related helpers. These functions return
    [`CompletelyReducibleBundle`](@ref) objects whose internal decomposition
    into `IrrepLevi` components happens automatically.

    In particular, `IrrepLevi` is **not** intended as the default user-facing
    constructor for bundles. When constructing a bundle directly from ambient
    weights, prefer [`CompletelyReducibleBundle`](@ref)`(X, coeffs)` or
    [`CompletelyReducibleBundle`](@ref)`(X, weights)` over manually creating
    `IrrepLevi` components first.

    The situations where a user encounters `IrrepLevi` directly are:
    - Inspecting [`components`](@ref)`(E)` to see the irreducible summands of
      a bundle
    - Advanced computations that need the central/semisimple coordinates

!!! note "Caching of tensor products"
    Tensor products of `IrrepLevi` decompose via Semisimple.jl character arithmetic,
    which involves computing tensor product branching rules from Weyl character
    data. This is the hottest path in the package and is **LRU-cached**: see
    [`PartialFlagVarieties.configure_caches!`](@ref) for tuning the cache budget.

## Type

```@docs
IrrepLevi
```

## Constructors

```@docs
IrrepLevi(::MarkedDynkinType, ::AbstractVector{<:Integer})
```

## Accessors

```@docs
central_part
semisimple_part
p_dominant_weight
degree(::IrrepLevi)
is_p_dominant
to_ambient_weight
```

## Monoidal operations

```@docs
tensor_product(::IrrepLevi, ::IrrepLevi)
dual(::IrrepLevi)
exterior_power(::IrrepLevi, ::Integer)
symmetric_power(::IrrepLevi, ::Integer)
```
