# CompletelyReducibleBundle

Completely reducible equivariant vector bundles on ``G/P``.

An **equivariant vector bundle** on ``G/P`` corresponds to a representation
of the parabolic subgroup ``P``. Because ``P`` is not reductive (it has a
nontrivial unipotent radical ``U``), such a representation need not be
completely reducible. However, the Levi factor ``L = P/U`` **is** reductive,
so every finite-dimensional ``L``-module is completely reducible — a direct sum
of irreducibles. Passing from ``P``-representations to ``L``-representations
amounts to taking the **semisimplification** of the bundle: one forgets the
extensions between the composition factors.

This package implements only semisimplified bundles:
a [`CompletelyReducibleBundle`](@ref) stores a formal sum of
[`IrrepLevi`](@ref) components. This is adequate for computing sheaf
cohomology (which depends only on the composition factors, not on extensions)
and tensor algebra, but does not capture filtration data. For the latter, see
[`FilteredBundle`](@ref).

!!! note "Design note"
    The abstract type [`Bundle`](@ref) has two concrete subtypes:
    `CompletelyReducibleBundle` (the semisimplification — a formal direct sum)
    and [`FilteredBundle`](@ref) (a bundle with a filtration by subbundles,
    retaining the ordering). Most user-facing operations produce
    `CompletelyReducibleBundle`; filtered bundles arise from
    [`filtered_tangent_bundle`](@ref) and its derived operations.

    Tensor products and exterior powers of `CompletelyReducibleBundle` are
    computed by decomposing into Lie.jl character arithmetic, then
    extracting dominant weights. Identical summands are **deduplicated**:
    components that appear with multiplicity are stored once with a
    multiplicity count, keeping the representation compact.

## Types

```@docs
Bundle
CompletelyReducibleBundle
```

## Constructors

```@docs
CompletelyReducibleBundle(::PartialFlagVariety, ::AbstractVector{<:Integer})
CompletelyReducibleBundle(::PartialFlagVariety, ::AbstractVector{<:AbstractVector{<:Integer}})
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
O
zero_bundle
line_bundle
tangent_bundle
cotangent_bundle
canonical_bundle
anticanonical_bundle
det
determinant
```

## Operations

```@docs
dual(::CompletelyReducibleBundle)
tensor_product(::CompletelyReducibleBundle, ::CompletelyReducibleBundle)
direct_sum
exterior_power(::CompletelyReducibleBundle, ::Integer)
symmetric_power(::CompletelyReducibleBundle, ::Integer)
twist
```

## Arithmetic

```@docs
Base.:*(::Integer, ::CompletelyReducibleBundle)
Base.iszero(::CompletelyReducibleBundle)
```
