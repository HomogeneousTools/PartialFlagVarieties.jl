# IrrepLevi

Irreducible representations of the Levi subgroup.

## Type

```@docs
IrrepLevi
```

## Accessors

```@docs
central_part
semisimple_part
fiber_dimension
to_ambient_weight
```

## Monoidal operations

```@docs
tensor_product(::IrrepLevi, ::IrrepLevi)
dual(::IrrepLevi)
exterior_power(::IrrepLevi, ::Int)
symmetric_power(::IrrepLevi, ::Int)
```
