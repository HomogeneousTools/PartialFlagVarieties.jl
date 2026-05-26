# Zero Loci

Zero loci of sections of equivariant bundles on ``G/P``.

Given an equivariant bundle ``E`` of rank ``r`` on ``X = G/P`` and a regular section
``s \in H^0(X, E)``, the **zero locus** ``Z = Z(s) \subset X`` is a smooth
subvariety of codimension ``r``. Its structure sheaf is resolved by the
**Koszul complex**

```math
0 \to \bigwedge^r E^* \to \cdots \to \bigwedge^2 E^* \to E^* \to \mathcal{O}_X \to \mathcal{O}_Z \to 0.
```

Twisting by any equivariant bundle ``F`` on ``X`` and taking the long exact
sequence in cohomology, one recovers ``H^*(Z, F|_Z)`` from computable data on
``X``. This is the basis for all invariants computed on zero loci:
cohomology of restrictions, Hodge numbers, Hilbert polynomials, and the
Calabi–Yau / Fano classification.

The package takes the geometric point of view that a **regular section has been
chosen**. The constructor `zero_locus(E)` therefore studies the geometry that
would result from a regular section of `E`; it does not prove existence of such
a section.

See the [mathematical background](../math.md#Koszul-resolution-and-zero-loci)
for details on the Koszul resolution and the
[conormal filtration](../math.md#The-conormal-filtration) that reduces Hodge
number computations to data on the ambient variety.

If you need compact serialization for datasets or classifications, see
[Labels](labels.md).

!!! warning "Lefschetz caveat for higher-rank bundles"
    The **Lefschetz hyperplane theorem** guarantees ``\mathrm{Pic}(X) \xrightarrow{\sim} \mathrm{Pic}(Z)``
    only when ``Z`` is a hypersurface (``\mathrm{rank}\, E = 1``). For zero loci
    of higher-rank bundles, the Picard rank of ``Z`` can exceed that of ``X``,
    so `hodge_numbers_symbolic` may leave ``h^{1,1}`` as an underdetermined
    symbolic variable — **this is correct behaviour**, not a bug.

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

## Koszul complex

```@docs
koszul_terms
```

## Cohomology on restrictions

```@docs
cohomology_on_restriction
cohomology_on_restriction_symbolic
euler_characteristic(::ZeroLocus, ::CompletelyReducibleBundle)
```

## Cohomological invariants

```@docs
hodge_numbers_symbolic
hodge_numbers_les
hilbert_polynomial(::ZeroLocus)
fano_index(::ZeroLocus)
euler_characteristic_tangent_bundle
```

## Classification predicates

```@docs
is_calabi_yau
is_strict_calabi_yau
is_strongly_fano
```
