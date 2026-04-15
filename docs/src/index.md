# PartialFlagVarieties.jl

A Julia package for computing with **partial flag varieties** ``G/P`` using
[Lie.jl](https://github.com/HomogeneousTools/Lie.jl).

## Features

- **Partial flag varieties** for all simple Lie types
  (``\mathrm{A}``–``\mathrm{G}``, including ``\mathrm{E}_6``,
  ``\mathrm{E}_7``, ``\mathrm{E}_8``, ``\mathrm{F}_4``, ``\mathrm{G}_2``)
- **Named constructors**: `Gr`, `OGr`, `SGr`, `LGr`, `IGr`,
  `projective_space`, `quadric`, `flag_variety`, `cayley_plane`,
  `freudenthal_variety`, `adjoint_variety`, `coadjoint_variety`
- **Equivariant vector bundles**: structure sheaf, tangent/cotangent,
  canonical, line bundles, exterior/symmetric powers, tensor products,
  duals, twists, determinants
- **Universal bundles**: tautological subbundles, quotient bundles, spinor
  bundles on quadrics
- **Filtered bundles**: tangent bundle filtration by root height, with
  induced filtrations on exterior/symmetric powers, duals, and tensor
  products
- **Sheaf cohomology** via the Borel–Weil–Bott theorem (character-valued
  and dimension-valued)
- **Hilbert polynomials** of equivariant bundles and zero loci
- **Zero loci** of sections of equivariant bundles: Koszul resolutions,
  restriction cohomology, Calabi–Yau detection
- **Hodge numbers**, **twisted Hodge numbers**, and **Hochschild
  cohomology** with polyvector parallelogram display
- **Symbolic Hodge computation** for zero loci, using long exact sequences,
  Serre duality cross-constraints, and Akizuki–Nakano vanishing
- **Exceptional collections**: Beilinson on ``\mathbb{P}^n``, Kapranov
  on quadrics, Kapranov–Orlov on Grassmannians, Schur functors
- **Topological invariants**: dimension, Euler characteristic, Betti
  numbers, Picard rank, Fano index, classification predicates
  (minuscule, cominuscule, adjoint, coadjoint)

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
``Z(L)^\circ``) and a semisimple part (highest weight of ``[L,L]``).

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
  "api/hochschild.md",
  "api/zero_loci.md",
  "api/koszul.md",
  "api/constructions.md",
  "api/exceptional_collections.md",
]
```
