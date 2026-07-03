# PartialFlagVarieties.jl

A Julia package for computing with **partial flag varieties** ``\mathrm{G}/\mathrm{P}`` using
[Semisimple.jl](https://github.com/HomogeneousTools/Semisimple.jl).

The package is designed for three closely related tasks:

1. **Construct homogeneous spaces ``\mathrm{G}/\mathrm{P}``.** Build partial flag varieties from
   Dynkin data or from named families such as Grassmannians, quadrics, and
   exceptional varieties.
2. **Compute with equivariant bundles.** Build semisimplified equivariant
   bundles and compute their sheaf cohomology via Borel–Weil–Bott.
3. **Study zero loci.** Form the zero locus of a regular section, run its
   Koszul resolution, and extract Hodge, Hochschild, Calabi–Yau, and Fano data.

## Start here

Different readers want different kinds of documentation:

| If you want to... | Start with... |
|:------------------|:--------------|
| get a first computation running | [Getting Started](getting_started.md) |
| understand the package conventions | [Conventions & Notation](conventions.md) |
| find a task-oriented example | [Common Workflows](workflows.md) |
| understand the mathematics behind the algorithms | [Mathematical Background](math.md) |
| look up a specific type or function | the pages under **API Reference** |

The API is intentionally layered. Most users should begin with:

- [`PartialFlagVariety`](@ref) and the [named constructors](api/constructions.md),
- [`CompletelyReducibleBundle`](@ref) together with `line_bundle`,
  `tangent_bundle`, `direct_sum`, `exterior_power`, and `twist`,
- [`cohomology`](@ref) / [`dimensions`](@ref),
- [`ZeroLocus`](@ref), [`hodge_numbers`](@ref), and [`hochschild_cohomology`](@ref).

Types such as [`IrrepLevi`](@ref) and the symbolic Koszul solvers are
documented because they matter for advanced use and for understanding the
internals, but they are not the recommended entry point.

## Quick taste

```julia
using PartialFlagVarieties

X = Gr(2, 5)
T = tangent_bundle(X)

dimension(X)              # 6
dimensions(T)[0]          # h^0(X, T_X)

Z = zero_locus(O(X, 1))
hodge_numbers(Z)
```

The most important indexing convention is that [`Cohomology`](@ref) objects are
**0-indexed**: `H[0]` means ``\mathrm{H}^0``.

## Mental model

At the user level, most computations follow this pattern:

1. Construct the ambient variety `X::PartialFlagVariety`.
2. Build a bundle `E::CompletelyReducibleBundle` on `X`.
3. Compute either:
   - cohomology on `X`,
   - or geometry/cohomology of `zero_locus(E)`.

Internally, each `PartialFlagVariety` wraps a [`MarkedDynkinType`](@ref) with
runtime Dynkin data and cached structural invariants. Equivariant bundles are
handled through the Levi factor, so the main bundle type stores the
**semisimplification** as a direct sum of irreducible Levi representations.

## Contents

```@contents
Pages = [
  "getting_started.md",
  "conventions.md",
  "workflows.md",
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
  "api/labels.md",
  "api/exceptional_collections.md",
  "api/cache_config.md",
]
```

## Authors

- Pieter Belmans <[pbelmans@uu.nl](mailto:pbelmans@uu.nl)>
- Javier Fernández Píriz <[javier.fernandezpiriz@uni.lu](mailto:javier.fernandezpiriz@uni.lu)>
