# MarkedDynkinType

The core type encoding a partial flag variety ``\mathrm{G}/\mathrm{P}_I`` and the
Levi factor of the parabolic subgroup ``\mathrm{P}_I``.

A `MarkedDynkinType` stores the ambient Dynkin type of ``\mathrm{G}`` together with a
set of **marked nodes** ``I \subseteq \{1, \ldots, r\}`` that specify the
parabolic ``\mathrm{P}_I``. Structural invariants — the Levi type, decomposition matrix,
root data, dimension, Euler characteristic, Betti numbers — are computed on
demand and cached per marked Dynkin type.

## Type

```@docs
MarkedDynkinType
```

## Accessors

```@docs
dynkin_type
marked_nodes
unmarked_nodes
central_rank
levi_rank
levi_type
is_borel
```

## Decomposition matrix

The **decomposition matrix** ``M`` performs the change of basis from the
fundamental weight coordinates of ``\mathrm{G}`` to the (center + semisimple) coordinate
system of the Levi factor ``\mathrm{L} = \operatorname{Z}(\mathrm{L})^\circ \cdot [\mathrm{L}, \mathrm{L}]``. Concretely, for
the marked nodes ``I`` the matrix ``M`` has ``|I|`` rows from the inverse
Cartan matrix ``C^{-1}`` and identity rows for the unmarked nodes, after
reordering via `levi_permutation`.

```@docs
decomposition_matrix
decomposition_matrix_inv
central_scaling_factor
levi_permutation
```

## Roots

```@docs
positive_nonparabolic_roots
positive_parabolic_roots
tangent_weights
```

## Display

```@docs
marked_dynkin_diagram
```

## Semisimple.jl utilities

Working with Dynkin types themselves is
[Semisimple.jl](https://github.com/HomogeneousTools/Semisimple.jl)'s job, and this package
re-exports the three functions used most often alongside a `MarkedDynkinType`:
`parse_dynkin_type` reads a type such as `"A2xB3"` from a string, while `sub_dynkin_type`
and `sub_dynkin_ordering` identify the sub-diagram induced on a set of nodes, which is what
[`levi_type`](@ref) and [`levi_permutation`](@ref) are built on.
