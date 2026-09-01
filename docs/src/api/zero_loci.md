# Zero Loci

Zero loci of sections of equivariant bundles on ``\mathrm{G}/\mathrm{P}``.

Given an equivariant bundle ``\mathcal{E}`` of rank ``r`` on ``X = \mathrm{G}/\mathrm{P}`` and a regular section
``s \in \mathrm{H}^0(X, \mathcal{E})``, the **zero locus** ``Z = Z(s) \subset X`` is a smooth
subvariety of codimension ``r``. Its structure sheaf is resolved by the
**Koszul complex**

```math
0 \to \bigwedge\nolimits^r \mathcal{E}^\vee \to \cdots \to \bigwedge\nolimits^2 \mathcal{E}^\vee \to \mathcal{E}^\vee \to \mathcal{O}_X \to \mathcal{O}_Z \to 0.
```

Twisting by any equivariant bundle ``\mathcal{F}`` on ``X`` and taking the long exact
sequence in cohomology, one recovers ``\mathrm{H}^\bullet(Z, \mathcal{F}|_Z)`` from computable data on
``X``. This is the basis for all invariants computed on zero loci:
cohomology of restrictions, Hodge numbers, Hilbert polynomials, and the
Calabi–Yau / Fano classification.

The constructor `zero_locus(E)` records the ambient variety and defining bundle,
but it does not choose or store a section. It represents the formal geometry of
a regular zero locus of `E`, assuming that a regular section exists.
Accordingly, equality compares the ambient variety and defining bundle. It does
not compare sections, test containment, or recognize abstractly isomorphic
subvarieties.

See the [mathematical background](../math.md#Koszul-resolution-and-zero-loci)
for details on the Koszul resolution and the
[conormal filtration](../math.md#The-conormal-filtration) that reduces Hodge
number computations to data on the ambient variety.

If you need compact serialization for datasets or classifications, see
[Labels](labels.md).

!!! warning "Lefschetz caveat for higher-rank bundles"
    The **Grothendieck–Lefschetz theorem** guarantees ``\mathrm{Pic}(X) \xrightarrow{\sim} \mathrm{Pic}(Z)``
    when ``Z`` is a smooth complete intersection of ample divisors of dimension
    ``\ge 3`` (in particular for ample hypersurfaces). For zero loci of
    higher-rank bundles that do **not** split as sums of ample line bundles,
    the Picard rank of ``Z`` can exceed that of ``X``, so `hodge_numbers` may
    leave ``\mathrm{h}^{1,1}`` as an underdetermined symbolic variable — **this is
    correct behaviour**, not a bug.

## Type

```@docs
ZeroLocus
zero_locus
```

## Accessors

```@docs
ambient_variety
defining_bundle
codimension
dimension(::ZeroLocus)
normal_bundle
conormal_bundle
```

## Relations

```@docs
is_sublocus
```

## Products

The product of formal zero loci is formed in the product ambient variety. Its
defining bundle is the [`external_direct_sum`](@ref) of the two defining
bundles.

```@docs
product(::ZeroLocus, ::ZeroLocus, ::ZeroLocus...)
```

## Bundles on a zero locus

Bundles on a zero locus use the same public bundle interface as bundles on the
ambient partial flag variety. In particular, `variety`, `rank`, `dual`,
`tensor_product`, `exterior_power`, `symmetric_power`, `cohomology`, and
`euler_characteristic` work without a separate zero-locus API.

```@docs
PartialFlagVarieties.ZeroLocusBundle
```

Use `restrict(Z, F)` to restrict an ambient bundle. Intrinsic tangent and
cotangent bundles are constructed directly from `Z`:

```julia
X = projective_space(4)
Z = zero_locus(line_bundle(X, 5))

L = restrict(Z, line_bundle(X, 1))
TZ = tangent_bundle(Z)

cohomology(TZ)
cohomology(tensor_product(TZ, L))
cohomology(exterior_power(TZ, 2))
euler_characteristic(TZ)
```

Bundles on different zero loci can be combined on the product using
[`external_tensor_product`](@ref) and [`external_direct_sum`](@ref). Both
operations lift the ambient presentations, so they also support intrinsic and
composite bundles such as tangent bundles and their tensor powers.

Restriction also composes along formal complete intersections. If `Z1` is cut
out by `E` and `Z2` is cut out by `E ⊕ E′`, then `restrict(Z2, F)`
restricts a bundle `F` on `Z1` to `Z2`. This is precisely the
[`is_sublocus(Z2, Z1)`](@ref is_sublocus) relation. It uses [`is_summand`](@ref)
on the defining bundles because individual sections are not part of the data
model.

Internally, these bundles retain the terms of bounded ambient presentations;
the presentation maps are implicit. Tensor products totalize the terms.
Exterior and symmetric powers use the derived graded-power formula, so they
also work for arbitrary composite presentations rather than only for the
tangent and cotangent sequences.

!!! note "Generic bundle cohomology versus the Hodge engine"
    `cohomology(exterior_power(cotangent_bundle(Z), p))` evaluates that one
    bundle from its presentation. `hodge_numbers(Z)` uses the same conormal
    complexes in a specialized batch computation: it caches terms shared by
    different values of `p`, replaces the top exterior power by adjunction,
    and combines rows using Serre duality, Lefschetz constraints, exact Euler
    characteristics, and Kodaira–Akizuki–Nakano vanishing. It can therefore
    be faster and can determine entries that an isolated bundle computation
    correctly leaves symbolic.

```@docs
restrict
structure_sheaf(::ZeroLocus)
zero_bundle(::ZeroLocus)
line_bundle(::ZeroLocus, ::Integer)
tangent_bundle(::ZeroLocus)
cotangent_bundle(::ZeroLocus)
canonical_bundle(::ZeroLocus)
anticanonical_bundle(::ZeroLocus)
rank(::PartialFlagVarieties.ZeroLocusBundle)
dual(::PartialFlagVarieties.ZeroLocusBundle)
tensor_product(::PartialFlagVarieties.ZeroLocusBundle, ::PartialFlagVarieties.ZeroLocusBundle)
exterior_power(::PartialFlagVarieties.ZeroLocusBundle, ::Integer)
symmetric_power(::PartialFlagVarieties.ZeroLocusBundle, ::Integer)
det(::PartialFlagVarieties.ZeroLocusBundle)
euler_characteristic(::PartialFlagVarieties.ZeroLocusBundle)
cohomology(::PartialFlagVarieties.ZeroLocusBundle)
Base.iszero(::PartialFlagVarieties.ZeroLocusBundle)
```

## Koszul complex

```@docs
koszul_terms
```

## Cohomological invariants

```@docs
hodge_numbers_symbolic
hodge_numbers_les
hilbert_polynomial(::ZeroLocus, ::CompletelyReducibleBundle)
fano_index(::ZeroLocus)
```

## Classification predicates

```@docs
is_calabi_yau
is_strict_calabi_yau
is_strongly_fano
```

## Internals

```@docs
PartialFlagVarieties._AmbientBundlePresentation
PartialFlagVarieties._derived_power
```
