# MarkedDynkinType

The core type encoding a partial flag variety ``G/P_I`` and the
Levi factor of the parabolic subgroup ``P_I``.

A `MarkedDynkinType` stores the ambient Dynkin type of ``G`` together with a
set of **marked nodes** ``I \subseteq \{1, \ldots, r\}`` that specify the
parabolic ``P_I``. Structural invariants — the Levi type, decomposition matrix,
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
fundamental weight coordinates of ``G`` to the (center + semisimple) coordinate
system of the Levi factor ``L = \operatorname{Z}(L)^\circ \cdot [L, L]``. Concretely, for
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

## Lie.jl utilities

Convenience functions for working with Dynkin types, built on top of
[Lie.jl](https://github.com/ulsmart/Lie.jl).

```@docs
parse_dynkin_type
cartan_type
cartan_type_with_ordering
```
