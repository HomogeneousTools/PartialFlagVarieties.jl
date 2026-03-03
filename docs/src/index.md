# PartialFlagVarieties.jl

A Julia package for computing with **partial flag varieties** ``G/P`` using
compile-time specialization via [Lie.jl](../Lie.jl).

## Features

- **Type-level encoding**: Dynkin types and marked nodes are type parameters,
  enabling ``@generated`` functions for zero-cost abstractions.
- **Levi decomposition**: Automatic identification of the Levi subgroup type
  from the marked Dynkin diagram.
- **Bundle calculus**: Irreducible Levi representations, completely reducible
  bundles, and monoidal operations (tensor product, exterior/symmetric powers,
  dual, twist).
- **Sheaf cohomology**: Borel–Weil–Bott theorem applied to equivariant bundles,
  returning either Weyl characters or dimensions.
- **Named constructors**: `Gr(k,n)`, `OGr(k,n)`, `SGr(k,n)`, `projective_space(n)`,
  `quadric(n)`, `cayley_plane()`, and more.
- **Topological invariants**: Dimension, Euler characteristic, Betti numbers,
  Picard rank, classification predicates.

## Quick start

```julia
using PartialFlagVarieties

# The Grassmannian Gr(2, 5)
V = Gr(2, 5)
dimension(V)            # 6
euler_characteristic(V) # 10
betti_numbers(V)        # [1, 1, 1, 1, 1, 1, 1]

# Sheaf cohomology on ℙ⁴
MDT = MarkedDynkinType{TypeA{4}, (1,)}
L = line_bundle(MDT, 1)
H = dimensions(cohomology(MDT, L))
H[0]                    # 5 = dim H⁰(ℙ⁴, 𝒪(1))

# The Cayley plane
V = cayley_plane()
dimension(V)            # 16
euler_characteristic(V) # 27
```

## Design philosophy

The package follows the same compile-time specialization pattern as Lie.jl:

1. Each partial flag variety ``G/P`` is encoded as a `MarkedDynkinType{DT, Marked}`
   where both `DT` (Dynkin type) and `Marked` (crossed-out nodes) are **type parameters**.

2. Heavy mathematical computations (Cartan matrices, Levi decomposition, Betti numbers,
   special matrices) are performed at **compile time** via `@generated` functions,
   yielding zero runtime overhead for repeated queries.

3. Bundle operations use the Levi decomposition: each equivariant bundle decomposes
   as a direct sum of irreducible Levi representations `IrrepLevi{MDT}`, each of
   which has a central part (character of the center) and a semisimple part
   (highest weight of the derived subgroup).

## Contents

```@contents
Pages = [
  "math.md",
  "api/marked_dynkin_type.md",
  "api/partial_flag_variety.md",
  "api/irrep_levi.md",
  "api/bundle.md",
  "api/cohomology.md",
  "api/constructions.md",
]
```
