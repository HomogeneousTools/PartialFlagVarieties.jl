# PartialFlagVarieties.jl

A Julia package for computing with **partial flag varieties** ``G/P`` using
[Lie.jl](https://github.com/HomogeneousTools/Lie.jl).

The three main capabilities are:

1. **Equivariant vector bundles and sheaf cohomology.**  Construct completely
   reducible equivariant bundles (line bundles, tangent/cotangent, exterior
   and symmetric powers, tensor products, …) and compute their sheaf cohomology
   via the Borel–Weil–Bott theorem.

2. **Hodge numbers and Hochschild cohomology.**  Compute the full Hodge diamond
   of homogeneous varieties and their zero loci, including symbolic Hodge
   computation when the Koszul long exact sequences leave degrees of freedom.

3. **Zero loci of sections.**  Given a section of an equivariant bundle, form
   the zero locus, run its Koszul complex, compute Euler characteristics and
   Hodge numbers, and detect Calabi–Yau or Fano geometry.

## Quick start

```julia
using PartialFlagVarieties

# The Grassmannian Gr(2, 5)
V = Gr(2, 5)
dimension(V)            # 6
euler_characteristic(V) # 10
betti_numbers(V)        # [1, 1, 1, 1, 1, 1, 1]

# Sheaf cohomology on ℙ⁴
X = projective_space(4)
L = line_bundle(X, 1)
H = dimensions(cohomology(L))
H[0]                    # 5 = dim H⁰(ℙ⁴, 𝒪(1))

# Zero loci
X = Gr(2, 5)
Z = zero_locus(line_bundle(X, 1))
dimension(Z)            # 5

# The Cayley plane
V = cayley_plane()
dimension(V)            # 16
euler_characteristic(V) # 27
```

## Design

Each partial flag variety ``G/P`` is encoded as a
`PartialFlagVariety{MDT}` wrapping a `MarkedDynkinType` that stores the
Lie type and marked nodes as **runtime values**. Derived structural
invariants (Cartan matrices, Levi decomposition, Betti numbers,
decomposition matrices) are **cached on demand**.

Bundle operations use the Levi decomposition: each equivariant bundle
decomposes as a direct sum of irreducible Levi representations
`IrrepLevi`, each with a central part (character of the center
``\operatorname{Z}(L)^\circ``) and a semisimple part (highest weight of ``[L,L]``).

## Contents

```@contents
Pages = [
  "math.md",
  "api/marked_dynkin_type.md",
  "api/partial_flag_variety.md",
  "api/irrep_levi.md",
  "api/bundle.md",
  "api/filtered_bundle.md",
  "api/universal_bundles.md",
  "api/cohomology.md",
  "api/hodge.md",
  "api/zero_loci.md",
  "api/koszul.md",
  "api/constructions.md",
  "api/exceptional_collections.md",
  "conventions.md",
]
```
