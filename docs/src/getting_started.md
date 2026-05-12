# Getting Started

This page is the fastest way to become productive with the package. For the
mathematical conventions behind the API, see [Conventions & Notation](conventions.md).

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

This is currently necessary because `Base62`, `ZeroLocus62`, and `Lie.jl` are
not available from the General registry.

## The three basic objects

Most computations revolve around three kinds of values:

| Object | Meaning | Typical constructors |
|:-------|:--------|:---------------------|
| `PartialFlagVariety` | the ambient homogeneous space ``G/P`` | `Gr`, `quadric`, `partial_flag_variety` |
| `CompletelyReducibleBundle` | a semisimplified equivariant bundle on ``G/P`` | `line_bundle`, `tangent_bundle`, `direct_sum`, `exterior_power` |
| `ZeroLocus` | the zero locus of a regular section of such a bundle | `zero_locus(E)` |

If you are new to the package, prefer the **named constructors** and the
high-level bundle API. Types such as [`IrrepLevi`](@ref) are documented mainly
so you can understand internal bundle decompositions and inspect
[`components`](@ref)`(E)`, but they are not the normal starting point and most
users will rarely construct them directly.

## 1. Construct the ambient variety

```julia
using PartialFlagVarieties

X = Gr(2, 5)
dimension(X)             # 6
picard_rank(X)           # 1
betti_numbers(X)         # [b_0, b_2, ..., b_12]
```

When no named constructor exists, use [`partial_flag_variety`](@ref):

```julia
X = partial_flag_variety("E6", [1, 6])
```

Here the marked nodes are the crossed-out nodes of the Dynkin diagram; see
[Conventions & Notation](conventions.md#Marked-nodes).

## 2. Build bundles from the public API

```julia
X = Gr(2, 5)
U = universal_subbundle(X)
Q = universal_quotient_bundle(X)

E = direct_sum(line_bundle(X, 1), exterior_power(Q, 2))
rank_bundle(E)
```

The package works primarily with **completely reducible** bundles. This is the
right model for Borel–Weil–Bott computations. Filtered objects are available
separately when the order of the filtration matters.

## 3. Compute cohomology

```julia
X = projective_space(4)
L = line_bundle(X, 1)

Hchar = cohomology(L)     # character-valued
Hdim = dimensions(Hchar)  # dimension-valued

Hdim[0]                   # H^0
Hdim[1]                   # H^1
```

`Cohomology` objects use **0-based indexing**: `H[0]` means ``H^0``.

## 4. Form zero loci

```julia
X = Gr(2, 5)
Z = zero_locus(line_bundle(X, 1))

dimension(Z)
euler_characteristic(Z)
hodge_numbers(Z)
```

The constructor `zero_locus(E)` assumes that a **regular section** of `E`
exists and studies the geometry of its smooth zero locus. It does not solve the
existence problem for sections.

## 5. Move on to the guides

Once the basic pattern feels familiar, the next useful pages are:

1. [Conventions & Notation](conventions.md) for coordinates, indexing, and the marked-node convention.
2. [Common Workflows](workflows.md) for task-oriented examples.
3. [Mathematical Background](math.md) for the theorem-level explanation of the algorithms.
