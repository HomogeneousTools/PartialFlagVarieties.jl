# PartialFlagVarieties.jl

[![CI](https://github.com/pbelmans/PartialFlagVarieties.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/pbelmans/PartialFlagVarieties.jl/actions/workflows/CI.yml)

A Julia package for computing with **partial flag varieties** $G/P$:
equivariant bundles, sheaf cohomology via the Borel–Weil–Bott theorem,
Hodge numbers, Hochschild cohomology, and more.

Part of the [HomogeneousTools](https://homogeneous.tools) project, building on [Lie.jl](https://github.com/HomogeneousTools/Lie.jl).

## Features

- **Partial flag varieties** for all simple Lie types (A–G, including E₆, E₇, E₈, F₄, G₂)
- Named constructors: `Gr`, `OGr`, `SGr`, `LGr`, `projective_space`, `quadric`, `flag_variety`, `cayley_plane`, `freudenthal_variety`, `adjoint_variety`, `coadjoint_variety`
- **Equivariant vector bundles**: structure sheaf, tangent/cotangent, canonical, line bundles, exterior/symmetric powers, tensor products, duals, twists
- **Universal bundles**: tautological subbundle, quotient bundle, spinor bundles on quadrics
- **Filtered bundles**: tangent bundle filtration by root height (Lemma 2.1, [arXiv:1606.04076](https://arxiv.org/abs/1606.04076))
- **Sheaf cohomology** via the Borel–Weil–Bott theorem (dimension-valued and character-valued)
- **Hodge numbers**, **twisted Hodge numbers**, **Hochschild cohomology** with polyvector parallelogram display
- Type-level encoding with `@generated` functions for compile-time constants
- Bourbaki/Oscar conventions throughout

## Quick start

```julia
using PartialFlagVarieties

# Grassmannian Gr(2, 5)
X = Gr(2, 5)
dimension(X)         # 6
euler_characteristic(X)  # 10

# Tangent bundle cohomology
T = tangent_bundle(X)
dimensions(T)        # H⁰ = 24

# Exterior powers
E = exterior_power(T, 2)
dimensions(E)        # H⁰ = 276

# Universal bundles
U = universal_subbundle(X)
Q = universal_quotient_bundle(X)
rank_bundle(U)       # 2
rank_bundle(Q)       # 3

# Hodge numbers
H = hodge_numbers(X)  # diagonal: h^{p,p} = b_{2p}

# Hochschild cohomology (polyvector parallelogram)
P = hochschild_cohomology(projective_space(2))
P[1, 0]  # h⁰(T) = 8 = dim Aut(ℙ²)
P[2, 0]  # h⁰(∧²T) = 10
```

## Named varieties

```julia
projective_space(4)          # ℙ⁴
quadric(5)                   # Q⁵
Gr(2, 5)                     # Grassmannian
OGr(2, 7)                    # Orthogonal Grassmannian
SGr(2, 6)                    # Symplectic Grassmannian
LGr(3, 6)                    # Lagrangian Grassmannian
flag_variety(5, (1, 3))      # Fl(1,3; 5)
cayley_plane()               # E₆/P₁
freudenthal_variety()        # E₇/P₇
adjoint_variety(TypeE{6})    # E₆/P₂
```

## Spinor bundles

```julia
X = quadric(5)               # Q⁵ = B₃/P₁
S = spinor_bundle(X)         # rank 4

X = quadric(4)               # Q⁴ = D₃/P₁
Sp = spinor_bundle(X, :plus)   # rank 2
Sm = spinor_bundle(X, :minus)  # rank 2
```

## Filtered tangent bundle

```julia
X = Gr(2, 5)
F = filtered_tangent_bundle(X)
n_filtration_steps(F)   # number of graded pieces
total_bundle(F)         # the underlying tangent bundle
graded_pieces(F)        # the associated graded (CompletelyReducibleBundle[])
```

## Examples

The `examples/` directory contains standalone scripts:

- **`HochschildAffine.jl`** — Polyvector fields on affine cones of $G/P$ for all types up to rank 8, with progress bars, output files, and cache management
- **`BottVanishing.jl`** — Verification of Bott vanishing failure for (co)adjoint varieties

Run with:
```sh
julia --project=. examples/HochschildAffine.jl
julia --project=. examples/BottVanishing.jl
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/HomogeneousTools/Lie.jl")
Pkg.add(url="https://github.com/pbelmans/PartialFlagVarieties.jl")
```

Requires Julia ≥ 1.9.

## Documentation

Full documentation at [homogeneous.tools](https://homogeneous.tools).

## Testing

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```
