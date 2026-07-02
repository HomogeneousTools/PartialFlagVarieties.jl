# ═══════════════════════════════════════════════════════════════════════════════
#  CacheConfig.jl — LRU cache configuration and public management API
#
#  All global caches in PartialFlagVarieties use LRU{K,V} from LRUCache.jl
#  with memory-based eviction (by = Base.summarysize). This file provides:
#
#   - Preferences-based persistent configuration (read in __init__)
#   - Runtime API: configure_caches!, clear_caches!, cache_info
#
#  Included after all cache-defining source files (MarkedDynkinType, IrrepLevi,
#  CompletelyReducibleBundle, ZeroLoci) so that the LRU constants exist.
# ═══════════════════════════════════════════════════════════════════════════════

# ─── Apply Preferences overrides ─────────────────────────────────────────────

"""
    _apply_cache_preferences!()

Read cache configuration from Preferences and resize caches accordingly.
Called from `__init__()`.

Supported preferences (set via `Preferences.set_preferences!`):
- `cache_budget_bytes::Int`: total cache budget in bytes
- `tensor_product_cache_fraction::Float64`: fraction for tensor product cache
- `bwb_pair_cache_fraction::Float64`: fraction for BWB pair cache
- `structural_cache_fraction::Float64`: fraction for structural caches
"""
function _apply_cache_preferences!()
  budget = @load_preference("cache_budget_bytes", nothing)
  budget === nothing && return nothing

  configure_caches!(;
    budget=Int(budget),
    tensor_frac=Float64(
      @load_preference("tensor_product_cache_fraction", _DEFAULT_TENSOR_FRAC)
    ),
    bwb_frac=Float64(@load_preference("bwb_pair_cache_fraction", _DEFAULT_BWB_FRAC)),
    structural_frac=Float64(
      @load_preference("structural_cache_fraction", _DEFAULT_STRUCTURAL_FRAC)
    ),
  )
end

# ─── Public API ───────────────────────────────────────────────────────────────

"""
    configure_caches!(; budget::Integer=_default_cache_budget(),
                       tensor_frac::Real=0.70,
                       bwb_frac::Real=0.20,
                       structural_frac::Real=0.10)

Resize all LRU caches according to the given total memory budget and per-cache
fractions. The budget is in bytes; fractions should sum to ≤ 1.0.

 The structural fraction is split among the marked-Dynkin, tangent-reps,
 cotangent-reps, and cached cotangent exterior powers (40%, 25%, 25%, 10%
 respectively).

# Examples
```julia
# Use 2 GiB total for caches
configure_caches!(budget = 2 * 1024^3)

# Use 500 MiB with more budget for tensor products
configure_caches!(budget = 500 * 1024^2, tensor_frac = 0.80, bwb_frac = 0.15)
```
"""
function configure_caches!(;
  budget::Integer=_default_cache_budget(),
  tensor_frac::Real=_DEFAULT_TENSOR_FRAC,
  bwb_frac::Real=_DEFAULT_BWB_FRAC,
  structural_frac::Real=_DEFAULT_STRUCTURAL_FRAC,
)
  b = Int(budget)
  tf = Float64(tensor_frac)
  bf = Float64(bwb_frac)
  sf = Float64(structural_frac)

  resize!(_TENSOR_PRODUCT_CACHE; maxsize=_cache_maxsize(b, tf))
  resize!(_BWB_PAIR_CACHE; maxsize=_cache_maxsize(b, bf))
  resize!(_marked_dynkin_cache; maxsize=_cache_maxsize(b, sf * 0.4))
  resize!(_tangent_reps_cache; maxsize=_cache_maxsize(b, sf * 0.25))
  resize!(_cotangent_reps_cache; maxsize=_cache_maxsize(b, sf * 0.25))
  resize!(_COTANGENT_POWER_CACHE; maxsize=_cache_maxsize(b, sf * 0.1))
  nothing
end

"""
    clear_caches!()

Empty all PartialFlagVarieties and Semisimple.jl LRU caches.
Resets hit/miss counters.
"""
function clear_caches!()
  empty!(_marked_dynkin_cache)
  empty!(_tangent_reps_cache)
  empty!(_cotangent_reps_cache)
  empty!(_COTANGENT_POWER_CACHE)
  empty!(_TENSOR_PRODUCT_CACHE)
  empty!(_BWB_PAIR_CACHE)
  Semisimple.clear_all_caches!()
  nothing
end

"""
    cache_info() -> NamedTuple

Return a summary of all LRU cache statistics: hits, misses, current size,
and max size for each cache.

# Examples
```julia
julia> info = cache_info();

julia> info.tensor_product
CacheInfo(; hits=42, misses=10, currentsize=1234, maxsize=100000000)
```
"""
function cache_info()
  (
    marked_dynkin=LRUCache.cache_info(_marked_dynkin_cache),
    tangent_reps=LRUCache.cache_info(_tangent_reps_cache),
    cotangent_reps=LRUCache.cache_info(_cotangent_reps_cache),
    cotangent_powers=LRUCache.cache_info(_COTANGENT_POWER_CACHE),
    tensor_product=LRUCache.cache_info(_TENSOR_PRODUCT_CACHE),
    bwb_pair=LRUCache.cache_info(_BWB_PAIR_CACHE),
  )
end
