# Cohomology

Sheaf cohomology via the Borel–Weil–Bott theorem.

The central computation: given an equivariant bundle ``E`` on ``X = G/P``,
decompose ``E`` into irreducible Levi summands, lift each to an ambient weight,
and apply the BWB theorem to determine which cohomology degree (if any) is
nonzero. The result is a [`Cohomology`](@ref) object indexed from ``0``
(i.e. `H[0]` is ``\mathrm{H}^0``).

See the [mathematical background](../math.md#The-Borel–Weil–Bott-theorem) for
the theorem statement.

## Type

```@docs
Cohomology
```

## The Borel–Weil–Bott algorithm

```@docs
borel_weil_bott
```

## Computation

```@docs
cohomology
dimensions
```

## Derived invariants

```@docs
chi
euler_characteristic
hilbert_polynomial
Base.iszero(::Cohomology{BigInt})
```
