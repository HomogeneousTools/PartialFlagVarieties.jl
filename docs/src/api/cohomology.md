# Cohomology

Sheaf cohomology via the Borel–Weil–Bott theorem.

The central computation: given an equivariant bundle ``\mathcal{E}`` on ``X = \mathrm{G}/\mathrm{P}``,
decompose ``\mathcal{E}`` into irreducible Levi summands, lift each to an ambient weight,
and apply the BWB theorem to determine which cohomology degree (if any) is
nonzero. The result is a [`Cohomology`](@ref) object indexed from ``0``
(i.e. `H[0]` is ``\mathrm{H}^0``).

See the [mathematical background](../math.md#The-Borel–Weil–Bott-theorem) for
the theorem statement.

!!! note "0-based indexing"
    `Cohomology` objects are one of the few user-facing types in the package
    that do not follow ordinary Julia 1-based indexing. This is deliberate:
    `H[0]`, `H[1]`, ... mirror the mathematical notation ``\mathrm{H}^0``, ``\mathrm{H}^1``, ...

## Type

```@docs
Cohomology
```

## The Borel–Weil–Bott algorithm

```@docs
borel_weil_bott
```

## Computation

For a [`CompletelyReducibleBundle`](@ref) on a partial flag variety,
[`cohomology(E)`](@ref) returns dimensions by default; pass `characters=true`
to keep the full ambient representation in each cohomology group.

Cohomology of a [`FilteredBundle`](@ref) or [`ZeroLocusBundle`](@ref) is always
dimension-valued. Both types retain implicit maps: the differentials in the
filtration spectral sequence or ambient presentation may not determine every
dimension. In that case the result has [`AffineExpr`](@ref) entries. For a
zero-locus bundle, character-valued cohomology is moreover not geometrically
defined in general because the section cutting out the locus need not be
invariant under the ambient group.

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
