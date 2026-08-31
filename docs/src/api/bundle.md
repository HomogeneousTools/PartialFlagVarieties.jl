# CompletelyReducibleBundle

Completely reducible equivariant vector bundles on ``\mathrm{G}/\mathrm{P}``.

An **equivariant vector bundle** on ``\mathrm{G}/\mathrm{P}`` corresponds to a representation
of the parabolic subgroup ``\mathrm{P}``. Because ``\mathrm{P}`` is not reductive (it has a
nontrivial unipotent radical ``\mathrm{U}``), such a representation need not be
completely reducible. However, the Levi factor ``\mathrm{L} = \mathrm{P}/\mathrm{U}`` **is** reductive,
so every finite-dimensional ``\mathrm{L}``-module is completely reducible — a direct sum
of irreducibles. Passing from ``\mathrm{P}``-representations to ``\mathrm{L}``-representations
amounts to taking the **semisimplification** of the bundle: one forgets the
extensions between the composition factors.

A [`CompletelyReducibleBundle`](@ref) stores a formal sum of
[`IrrepLevi`](@ref) components. This is enough for Borel–Weil–Bott calculations
on each graded piece and for tensor algebra, but individual cohomology groups
of a nonsplit extension can depend on its connecting maps. For ordered
filtration data, see [`FilteredBundle`](@ref); bundles on formal zero loci use
ambient presentations as described under [Zero Loci](zero_loci.md).

!!! note "Design note"
    The abstract type [`Bundle`](@ref) has three concrete subtypes:
    `CompletelyReducibleBundle` (the semisimplification — a formal direct sum)
    [`FilteredBundle`](@ref) (a bundle with a filtration by subbundles,
    retaining the ordering), and the internal `ZeroLocusBundle` returned by
    zero-locus bundle constructors. Most ambient operations produce
    `CompletelyReducibleBundle`; filtered bundles arise from
    [`filtered_tangent_bundle`](@ref) and its derived operations.

    Tensor products and exterior powers of `CompletelyReducibleBundle` are
    computed by decomposing into Semisimple.jl character arithmetic, then
    extracting dominant weights. Repeated summands remain represented with
    their multiplicities; expensive power calculations group equal summands
    internally before expanding them.

## Types

```@docs
Bundle
CompletelyReducibleBundle
```

## Constructors

When constructing a bundle manually from ambient highest weights, prefer these
constructors over assembling `IrrepLevi` summands by hand.

```@docs
CompletelyReducibleBundle(::PartialFlagVariety, ::WeightLatticeElem)
CompletelyReducibleBundle(::PartialFlagVariety, ::AbstractVector{<:WeightLatticeElem})
CompletelyReducibleBundle(::PartialFlagVariety, ::AbstractVector{<:Integer})
CompletelyReducibleBundle(::PartialFlagVariety, ::AbstractVector{<:AbstractVector{<:Integer}})
E
```

## Accessors

```@docs
variety
components
n_components
rank(::CompletelyReducibleBundle)
```

## Positivity

```@docs
picard_degrees
is_ample_line_bundle
```

## Standard bundles

```@docs
structure_sheaf
O
zero_bundle
line_bundle
tangent_bundle
T
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
is_summand
exterior_power(::CompletelyReducibleBundle, ::Integer)
symmetric_power(::CompletelyReducibleBundle, ::Integer)
twist
```

## External operations

For bundles on different partial flag varieties, external operations pull both
bundles back to the product before applying the ordinary operation.

```@docs
external_tensor_product
external_direct_sum
```

## Arithmetic

```@docs
Base.:*(::Integer, ::CompletelyReducibleBundle)
Base.iszero(::CompletelyReducibleBundle)
```
