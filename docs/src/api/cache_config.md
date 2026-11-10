# Cache Configuration

Runtime management of the LRU caches that accelerate repeated computations.

The package caches five categories of data using in-memory LRU eviction
(via [LRUCache.jl](https://github.com/JuliaCollections/LRUCache.jl)):

| Cache | Contents |
|:------|:---------|
| `marked_dynkin` | Structural invariants of `MarkedDynkinType` values |
| `tangent_reps` | Tangent bundle decompositions into `IrrepLevi` |
| `cotangent_reps` | Cotangent bundle decompositions |
| `tensor_product` | Tensor products of `IrrepLevi` (the hot path) |
| `bwb_pair` | Borel–Weil–Bott pair lookups |

By default the total cache budget is determined automatically from available
system memory. Use [`PartialFlagVarieties.configure_caches!`](@ref) to override, or set persistent
preferences via `Preferences.set_preferences!`.

## API

```@docs
PartialFlagVarieties.configure_caches!
PartialFlagVarieties.clear_caches!
PartialFlagVarieties.cache_info
```

Projective-duality computations share a cache for the chosen embedding,
the resulting invariants, and the section Euler characteristics used to find them.
Its key includes the line bundle, so changing a stored polarization selects
the appropriate computation. It is included in `configure_caches!`,
`clear_caches!`, and `cache_info().projective_dual`.
