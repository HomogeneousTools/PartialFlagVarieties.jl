# PartialFlagVarieties.jl

[![Tests](https://github.com/HomogeneousTools/PartialFlagVarieties.jl/actions/workflows/test.yml/badge.svg)](https://github.com/HomogeneousTools/PartialFlagVarieties.jl/actions/workflows/test.yml)
[![Docs](https://img.shields.io/badge/docs-homogeneous.tools/PartialFlagVarieties.jl-blue)](https://homogeneous.tools/PartialFlagVarieties.jl)
[![Release](https://img.shields.io/github/v/release/HomogeneousTools/PartialFlagVarieties.jl?color=green)](https://github.com/HomogeneousTools/PartialFlagVarieties.jl/releases)

A Julia package for computing with **partial flag varieties** $G/P$:
equivariant vector bundles, sheaf cohomology via the Borel–Weil–Bott theorem,
zero loci, Hodge numbers, Hochschild cohomology, exceptional collections,
and more.

Part of the [HomogeneousTools](https://homogeneous.tools) project, building on [Lie.jl](https://github.com/HomogeneousTools/Lie.jl).

## Features

- **Partial flag varieties** for all simple Lie types (A–G, including $\mathrm{E}_6$, $\mathrm{E}_7$, $\mathrm{E}_8$, $\mathrm{F}_4$, $\mathrm{G}_2$)
- Named constructors: `Gr`, `OGr`, `SGr`, `LGr`, `IGr`, `projective_space`, `quadric`, `flag_variety`, `cayley_plane`, `freudenthal_variety`, `adjoint_variety`, `coadjoint_variety`
- **Equivariant vector bundles**: structure sheaf, tangent/cotangent, canonical, line bundles, exterior/symmetric powers, tensor products, duals, twists, determinants
- **Universal bundles**: tautological subbundles, quotient bundles, spinor bundles on quadrics
- **Filtered bundles**: tangent bundle filtration by root height, with induced filtrations on exterior/symmetric powers, duals, and tensor products
- **Sheaf cohomology** via the Borel–Weil–Bott theorem (both dimension-valued and character-valued)
- **Hilbert polynomials** of equivariant bundles and zero loci
- **Zero loci** of sections of equivariant bundles, with Koszul resolutions, restriction cohomology, and Calabi–Yau detection
- **Hodge numbers**, **twisted Hodge numbers**, **Hochschild cohomology** with polyvector parallelogram display
- **Symbolic Hodge computation** for zero loci using long exact sequences, Serre duality cross-constraints, and Akizuki–Nakano vanishing
- **Exceptional collections**: Beilinson on $\mathbb{P}^n$, Kapranov on quadrics, Kapranov–Orlov on Grassmannians, Schur functors
- Topological invariants: dimension, Euler characteristic, Betti numbers, Picard rank, Fano index, classification predicates (minuscule, cominuscule, adjoint, coadjoint)
- Runtime marked Dynkin data with cached structural invariants
- Bourbaki/Oscar conventions throughout

## Quick start

```julia
using PartialFlagVarieties

# Grassmannian Gr(2, 5)
X = Gr(2, 5)
dimension(X)             # 6
euler_characteristic(X)  # 10

# Tangent bundle cohomology
T = tangent_bundle(X)
dimensions(T)            # H⁰ = 24

# Exterior powers
E = exterior_power(T, 2)
dimensions(E)            # H⁰ = 276

# Universal bundles
U = universal_subbundle(X)
Q = universal_quotient_bundle(X)
rank_bundle(U)           # 2
rank_bundle(Q)           # 3

# Hodge numbers
H = hodge_numbers(X)     # diagonal: h^{p,p} = b_{2p}

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

## Zero loci

```julia
X = Gr(2, 5)
L = line_bundle(X, 1)
Z = zero_locus(L)            # hypersurface in Gr(2,5)
dimension(Z)                 # 5
is_calabi_yau(Z)             # false

# Hodge numbers of zero loci
hodge_numbers(Z)
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

# Exterior/symmetric powers inherit the filtration
F2 = exterior_power(F, 2)   # ∧² with induced filtration
S2 = symmetric_power(F, 2)  # Sym² with induced filtration
```

## Exceptional collections

```julia
X = projective_space(3)
Es = beilinson_collection(X)               # ⟨𝒪, 𝒪(1), 𝒪(2), 𝒪(3)⟩
is_full_exceptional_sequence(Es, X)        # true

X = quadric(4)
Es = kapranov_collection(X)               # spinor bundles + line bundles

X = Gr(2, 4)
Es = kapranov_bundles_grassmannian(X)      # Schur functors of U∨
is_strong_exceptional_sequence(Es)         # true
```

## Examples

The `examples/` directory contains standalone scripts:

- **`HochschildAffine.jl`** — Euler characteristics of polyvector fields $\chi(\bigwedge^p T_{G/P})$ for all types up to rank 8
- **`BottVanishing.jl`** — Verification of Bott vanishing failure for (co)adjoint varieties
- **`CICY3-1606.04076.jl`** — Calabi–Yau threefold classification in exceptional flag varieties
- **`CICY3-1607.07821.jl`** — Complete-intersection Calabi–Yau threefolds on Grassmannians
- **`Kuechle.jl`** — Küchle's classification of Fano fourfolds in Grassmannians
- **`FanoThreefolds.jl`** — Hodge diamonds and polyvector parallelograms for homogeneous Fano threefolds
- **`FanoFourfolds.jl`** — Hodge number computation for Fano fourfolds from JSON dataset
- **`ExceptionalCollections.jl`** — Exceptional collection verification on various varieties
- **`Hyperkaehler.jl`** — Hodge numbers of hyperkähler fourfolds of $\mathrm{K3}^{[2]}$-type
- **`LinearSections.jl`** — Linear sections of Grassmannians and homological projective duality

Run with:
```sh
julia --project=. examples/HochschildAffine.jl
julia --project=. examples/BottVanishing.jl
```

## Installation

```julia
using Pkg
Pkg.add([
  Pkg.PackageSpec(url="https://github.com/HomogeneousTools/Base62"),
  Pkg.PackageSpec(url="https://github.com/HomogeneousTools/ZeroLocus62", subdir="julia"),
  Pkg.PackageSpec(url="https://github.com/HomogeneousTools/Lie.jl"),
  Pkg.PackageSpec(url="https://github.com/HomogeneousTools/PartialFlagVarieties.jl"),
])
```

`Base62`, `ZeroLocus62`, and `Lie.jl` are not currently available from the
General registry, so they must be added explicitly in a clean environment.
Listing them in `Project.toml` records the dependency graph, but it does not
teach `Pkg.add(url=...)` where to fetch unregistered dependencies when
`PartialFlagVarieties.jl` is installed into another environment.

Requires Julia ≥ 1.10.

## Documentation

Full documentation at [homogeneous.tools](https://homogeneous.tools).

## Testing

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Doctests

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate(); using Documenter, PartialFlagVarieties; doctest(PartialFlagVarieties)'
```

## Formatting

This repository uses [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl)
with the Blue style configured in `.JuliaFormatter.toml`.

```sh
julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("JuliaFormatter"); using JuliaFormatter; format(".")'
```
