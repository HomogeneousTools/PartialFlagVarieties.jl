# Zero Loci

Zero loci of sections of equivariant bundles on ``G/P``.

Given an equivariant bundle ``E`` on ``X = G/P`` and a regular section
``s \in H^0(X, E)``, the **zero locus** ``Z = Z(s) \subset X`` has
codimension equal to ``\operatorname{rank} E``. The Koszul complex

```math
0 \to \bigwedge^r E^* \to \cdots \to E^* \to \mathcal{O}_X \to \mathcal{O}_Z \to 0
```

resolves ``\mathcal{O}_Z`` and is the basis for computing cohomology on ``Z``.

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
```

## Hodge numbers

```@docs
hodge_numbers_symbolic
hodge_numbers_les
```

## Classification predicates

```@docs
is_calabi_yau
is_calabi_yau_candidate
is_weak_fano
```
