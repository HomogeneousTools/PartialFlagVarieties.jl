# ═══════════════════════════════════════════════════════════════════════════════
#  ZeroLoci.jl — Zero loci of sections of equivariant bundles
#
#  A ZeroLocus represents the zero locus Z(s) of a regular section
#  s ∈ H⁰(X, E) of an equivariant bundle E on X = G/P.
#
#  Provides:
#   - Euler characteristic of restrictions (always exact)
#   - Cohomology of restrictions via Koszul + LES
#   - Hodge numbers (when determined by the LES)
#   - Calabi–Yau detection
# ═══════════════════════════════════════════════════════════════════════════════

export ZeroLocus
export zero_locus, ambient_variety, defining_bundle
export factors, n_factors
export codimension, normal_bundle, conormal_bundle
export koszul_terms, cohomology_on_restriction, cohomology_on_restriction_symbolic
export is_calabi_yau, is_strict_calabi_yau
export is_strongly_fano
export fano_index
export hilbert_polynomial
export hodge_numbers_symbolic
export hodge_numbers_les
export euler_characteristic_tangent_bundle
export tangent_cohomology

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ZeroLocus

The zero locus ``Z(s)`` of a regular section ``s \\in \\mathrm{H}^0(X, \\mathcal{E})`` of an
equivariant bundle ``\\mathcal{E}`` on the partial flag variety ``X = \\mathrm{G}/\\mathrm{P}``.

Assumes the section is regular, so ``\\dim Z = \\dim X - \\mathrm{rank}(\\mathcal{E})``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> dimension(Z)
3

julia> euler_characteristic(Z)   # topological χ of the quintic Calabi–Yau
-200
```
"""
mutable struct ZeroLocus
  const ambient::PartialFlagVariety
  const defining_bundle::CompletelyReducibleBundle
  # Exterior powers of the dual bundle: koszul_wedges[i+1] = ∧ⁱE*.
  # Populated lazily on the first call that needs them; reused thereafter.
  koszul_wedges::Union{Nothing,Vector{CompletelyReducibleBundle}}
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Constructors
# ═══════════════════════════════════════════════════════════════════════════════

"""
    zero_locus(E::CompletelyReducibleBundle) -> ZeroLocus

Construct the zero locus of a regular section of the equivariant
bundle ``\\mathcal{E}``.  Requires ``\\mathrm{rank}(\\mathcal{E}) \\le \\dim(X)``.

This constructor assumes such a regular section exists; it does not try to
prove existence or regularity.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> dimension(Z)
3
```
"""
function zero_locus(E::CompletelyReducibleBundle)
  X = E.variety
  r = Int(rank(E))
  d = dimension(X)
  r <= d || throw(ArgumentError(
    "Bundle rank $r exceeds ambient dimension $d."
  ))
  ZeroLocus(X, E, nothing)
end

"""
    _product_zero_locus(Z1, Z2) -> ZeroLocus

Construct the product of two zero loci by taking the product ambient and the
direct sum of the defining bundles lifted from the two factors.
"""
function _product_zero_locus(Z1::ZeroLocus, Z2::ZeroLocus)
  X1 = ambient_variety(Z1)
  X2 = ambient_variety(Z2)
  X = product(X1, X2)
  E1 = _lift_bundle_to_product(X, defining_bundle(Z1), 0)
  E2 = _lift_bundle_to_product(X, defining_bundle(Z2), rank(dynkin_type(X1)))
  zero_locus(direct_sum(E1, E2))
end

"""
    product(Z1::ZeroLocus, Z2::ZeroLocus, Zs::ZeroLocus...) -> ZeroLocus

Construct the product of zero loci.

If `Z_i ⊂ X_i` is cut out by a regular section of `E_i`, then the product is
realized as the zero locus in `X_1 × X_2 × ...` of the direct sum of the
lifted bundles pulled back from each factor. This is also available through the
`*` operator.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> Z = product(zero_locus(line_bundle(projective_space(1), 1)),
                   zero_locus(line_bundle(projective_space(2), 1)));

julia> dimension(Z)
1
```
"""
function product(Z1::ZeroLocus, Z2::ZeroLocus, Zs::ZeroLocus...)
  Z = _product_zero_locus(Z1, Z2)
  for W in Zs
    Z = _product_zero_locus(Z, W)
  end
  Z
end

Base.:*(Z1::ZeroLocus, Z2::ZeroLocus) = product(Z1, Z2)

# ═══════════════════════════════════════════════════════════════════════════════
#  Product decomposition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    factors(Z::ZeroLocus) -> Vector{ZeroLocus}

Decompose `Z` into its product factors. Two ambient factors are *linked* when
some summand of the defining bundle is nontrivial on both; the connected
components of that relation give `Z = Z_1 × ⋯ × Z_k`, and this returns the
`Z_i` on their sub-product ambients. An ambient factor touched by no summand
splits off as the whole flag variety. Returns the singleton `[Z]` when `Z` is
irreducible in this "partition the ambient factors" sense.

A factor is kept when it is positive-dimensional or a set of `m ≥ 2` points (the
latter disconnects `Z` into `m` copies); a single reduced point is a Künneth
identity and is dropped. `n_factors` counts the kept factors.

`hodge_numbers`, `hochschild_cohomology` and `tangent_cohomology` recombine the
factors by the Künneth formula, which determines diamonds/parallelograms the
monolithic long-exact-sequence solver leaves symbolic. (The remaining
invariants — `euler_characteristic`, `hilbert_polynomial`, the anticanonical
degree — are already exact for a product via Koszul, so need no special path.)
"""
function factors(Z::ZeroLocus)
  ambient_factors = _mdt_to_factors(marked_dynkin_type(Z.ambient))
  length(ambient_factors) == 1 && return [Z]
  factor_ranks = Int[factor.rank for factor in ambient_factors]

  # Per summand, the per-ambient-factor weight blocks and the set of ambient
  # factors it actually twists (its support).
  summand_rows = [
    _weight_to_summand_row(p_dominant_weight(summand), factor_ranks)
    for summand in components(Z.defining_bundle)
  ]
  supports = [findall(block -> any(!iszero, block), row) for row in summand_rows]

  blocks = _connected_ambient_factors(supports, length(ambient_factors))
  length(blocks) == 1 && return [Z]

  parts = ZeroLocus[]
  for block in blocks
    subvariety = PartialFlagVariety(_factors_to_mdt(ambient_factors[block]))
    subtype = dynkin_type(subvariety)
    weights = [
      _summand_row_to_weight(row[block], subtype)
      for
      (row, support) in zip(summand_rows, supports) if support ⊆ block && !isempty(support)
    ]
    subbundle = if isempty(weights)
      zero_bundle(subvariety)
    else
      CompletelyReducibleBundle(subvariety, weights)
    end
    # An over-cut component makes the whole product empty; hand that back to the
    # monolithic solver rather than build an invalid (rank > dim) zero locus.
    rank(subbundle) > dimension(subvariety) && return [Z]
    push!(parts, zero_locus(subbundle))
  end
  # Keep positive-dimensional factors and multi-point (m ≥ 2) sets; drop single
  # reduced points, which are Künneth identities. `dimension >= 1` short-circuits,
  # so χ(𝒪) is evaluated only on 0-dimensional parts, where it equals h⁰·⁰, the
  # number of points (a cheap BWB alternating sum).
  filter(part -> dimension(part) >= 1 || euler_characteristic(part) >= 2, parts)
end

# Partition the ambient factors `1:n` into connected blocks, joining two factors
# whenever some bundle summand is supported on both (union–find with path
# compression). Untouched factors form singleton blocks.
function _connected_ambient_factors(supports, n)
  parent = collect(1:n)
  root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
  for support in supports, factor in support
    parent[root(factor)] = root(first(support))
  end
  blocks = [Int[] for _ in 1:n]
  for factor in 1:n
    push!(blocks[root(factor)], factor)
  end
  return filter(!isempty, blocks)
end

"""
    n_factors(Z::ZeroLocus) -> Int

Number of product factors of `Z` (see [`factors`](@ref)); equals `1` for an
irreducible `Z`. Invariant computations decompose via Künneth exactly when
`n_factors(Z) >= 2`.
"""
n_factors(Z::ZeroLocus) = length(factors(Z))

# ═══════════════════════════════════════════════════════════════════════════════
#  Accessors
# ═══════════════════════════════════════════════════════════════════════════════

"""
Return (and lazily compute) the cached vector `[∧⁰E*, ∧¹E*, …, ∧ʳE*]`.
All Koszul computations go through this accessor so the exterior powers
are computed at most once per `ZeroLocus`, regardless of how many times
different twists are requested.
"""
function _koszul_wedges!(Z::ZeroLocus)
  if Z.koszul_wedges === nothing
    r = Int(rank(Z.defining_bundle))
    E_dual = dual(Z.defining_bundle)
    Z.koszul_wedges = CompletelyReducibleBundle[exterior_power(E_dual, i) for i in 0:r]
  end
  Z.koszul_wedges
end

"""
    ambient_variety(Z::ZeroLocus) -> PartialFlagVariety

Return the ambient variety ``X`` containing the zero locus.
"""
ambient_variety(Z::ZeroLocus) = Z.ambient

"""
    defining_bundle(Z::ZeroLocus) -> CompletelyReducibleBundle

Return the bundle ``\\mathcal{E}`` whose section defines the zero locus.
"""
defining_bundle(Z::ZeroLocus) = Z.defining_bundle

"""
    codimension(Z::ZeroLocus) -> Int

Return the codimension of ``Z`` in ``X``, equal to ``\\mathrm{rank}(\\mathcal{E})``.
"""
codimension(Z::ZeroLocus)::Int = rank(Z.defining_bundle)

"""
    dimension(Z::ZeroLocus) -> Int

Return the dimension ``\\dim Z = \\dim X - \\mathrm{rank}(\\mathcal{E})`` of the zero locus.
"""
dimension(Z::ZeroLocus)::Int = dimension(Z.ambient) - codimension(Z)

# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul complex construction
# ═══════════════════════════════════════════════════════════════════════════════

"""
    koszul_terms(Z::ZeroLocus, F::CompletelyReducibleBundle)
      -> Vector{CompletelyReducibleBundle}

Return the terms of the twisted Koszul complex:
``[\\mathcal{F} \\otimes \\wedge^0 \\mathcal{E}^\\vee, \\mathcal{F} \\otimes \\wedge^1 \\mathcal{E}^\\vee, \\ldots,
   \\mathcal{F} \\otimes \\wedge^r \\mathcal{E}^\\vee]``
where ``\\mathcal{E}`` is the defining bundle and ``r = \\mathrm{rank}(\\mathcal{E})``.
"""
function koszul_terms(Z::ZeroLocus, F::CompletelyReducibleBundle)
  marked_dynkin_type(variety(F)) == marked_dynkin_type(Z.ambient) || throw(
    ArgumentError(
      "koszul_terms requires a bundle on the ambient variety of the zero locus."
    ),
  )
  CompletelyReducibleBundle[tensor_product(F, w) for w in _koszul_wedges!(Z)]
end

"""
    koszul_terms(Z::ZeroLocus) -> Vector{CompletelyReducibleBundle}

Return the terms of the (untwisted) Koszul complex for ``\\mathcal{O}_Z``.
"""
function koszul_terms(Z::ZeroLocus)
  koszul_terms(Z, structure_sheaf(Z.ambient))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Memory-efficient Koszul cohomology (bypass CRB construction)
# ═══════════════════════════════════════════════════════════════════════════════

# Cache: (a, b) → [(degree, dimension), ...] from BWB applied to tensor_product(a, b).
# Populated lazily; avoids repeated borel_weil_bott + degree calls.
const _BWB_PAIR_CACHE = let b = _default_cache_budget()
  LRU{Tuple{IrrepLevi,IrrepLevi},Vector{Pair{Int,BigInt}}}(;
    maxsize=_cache_maxsize(b, _DEFAULT_BWB_FRAC * 0.8),
    by=Base.summarysize,
  )
end
const _COTANGENT_POWER_CACHE = let b = _default_cache_budget()
  LRU{Tuple{MarkedDynkinType,Int},CompletelyReducibleBundle}(;
    maxsize=_cache_maxsize(b, _DEFAULT_STRUCTURAL_FRAC * 0.1),
    by=Base.summarysize,
  )
end

function _cotangent_power(X::PartialFlagVariety, j::Int)
  mdt = marked_dynkin_type(X)
  get!(_COTANGENT_POWER_CACHE, (mdt, j)) do
    exterior_power(cotangent_bundle(X), j)
  end
end

const _FILTERED_COTANGENT_POWER_CACHE = let b = _default_cache_budget()
  LRU{Tuple{MarkedDynkinType,Int},FilteredBundle}(;
    maxsize=_cache_maxsize(b, _DEFAULT_STRUCTURAL_FRAC * 0.1),
    by=Base.summarysize,
  )
end

"""
The ``j``-th exterior power of the filtered cotangent bundle, cached per
marked Dynkin type: the induced-filtration plethysm is the expensive part
of every filtered-ambient Hodge computation and only depends on ``(X, j)``.
"""
function _filtered_cotangent_power(X::PartialFlagVariety, j::Int)
  mdt = marked_dynkin_type(X)
  get!(_FILTERED_COTANGENT_POWER_CACHE, (mdt, j)) do
    exterior_power(dual(filtered_tangent_bundle(X)), j)
  end
end

"""
Compute and cache the BWB contributions `[(deg, dim), ...]` for
`tensor_product(a, b)`.  Returns an empty vector when all components are
acyclic.
"""
function _bwb_pair(a::IrrepLevi, b::IrrepLevi)
  get!(_BWB_PAIR_CACHE, _unordered_pair(a, b)) do
    by_degree = Dict{Int,BigInt}()
    for (c, mult) in _tensor_product_terms(a, b)
      λ = p_dominant_weight(c)
      bwb = borel_weil_bott(λ)
      bwb === nothing && continue
      deg, μ = bwb
      by_degree[deg] = get(by_degree, deg, BigInt(0)) + BigInt(mult) * BigInt(degree(μ))
    end
    Pair{Int,BigInt}[deg => dim for (deg, dim) in by_degree]
  end
end

"""
    _koszul_dimensions(Z::ZeroLocus, F::CompletelyReducibleBundle)
      -> Vector{Cohomology{BigInt}}

Compute dimension-valued cohomology for each Koszul term ``\\mathcal{F} ⊗ ∧^i \\mathcal{E}^\\vee``
without materialising intermediate `CompletelyReducibleBundle` objects.

Deduplicates F and wedge components into multiplicity dicts, then
computes tensor-product weight counts and applies BWB directly.
This avoids allocating ``O(N^2)`` `IrrepLevi` vectors that are
immediately re-deduplicated by `dimensions()`.
"""
function _koszul_dimensions(Z::ZeroLocus, F::CompletelyReducibleBundle)
  marked_dynkin_type(variety(F)) == marked_dynkin_type(Z.ambient) || throw(
    ArgumentError(
      "_koszul_dimensions requires a bundle on the ambient variety of the zero locus."
    ),
  )
  _koszul_dimensions(Z, _to_counts(F))
end

function _koszul_dimensions(Z::ZeroLocus)
  _koszul_dimensions(Z, structure_sheaf(Z.ambient))
end

"""
    _koszul_dimensions(Z::ZeroLocus, f_counts::Dict{IrrepLevi,Int})
      -> Vector{Cohomology{BigInt}}

Dict-accepting overload: skips CRB construction entirely.
"""
function _koszul_dimensions(
  Z::ZeroLocus, f_counts::Dict{IrrepLevi,Int}
)
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  _koszul_dimensions(Z, f_counts, wedge_counts)
end

function _koszul_dimensions(
  Z::ZeroLocus, f_counts::Dict{IrrepLevi,Int},
  wedge_counts::Vector{Dict{IrrepLevi,Int}},
)
  d = dimension(Z.ambient)

  result = Vector{Cohomology{BigInt}}(undef, length(wedge_counts))

  for (wi, w_counts) in enumerate(wedge_counts)
    entries = zeros(BigInt, d + 1)
    for (a, ma) in f_counts
      for (b, mb) in w_counts
        total = ma * mb
        for (deg, dim) in _bwb_pair(a, b)
          if 0 <= deg <= d
            entries[deg + 1] += total * dim
          end
        end
      end
    end
    result[wi] = Cohomology{BigInt}(entries, d)
  end

  result
end

"""
Dual of a multiplicity dict: map each IrrepLevi to its dual.
"""
function _dual_counts(f_counts::Dict{IrrepLevi,Int})
  result = Dict{IrrepLevi,Int}()
  for (component, mult) in f_counts
    dual_component = dual(component)
    result[dual_component] = get(result, dual_component, 0) + mult
  end
  result
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Euler characteristic (always exact)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    euler_characteristic(Z::ZeroLocus, F::CompletelyReducibleBundle) -> BigInt

Compute ``\\chi(Z, \\mathcal{F}|_Z) = \\sum_{i=0}^{r} (-1)^i \\chi(X, \\mathcal{F} \\otimes \\wedge^i \\mathcal{E}^\\vee)``.
This is always exact — no long exact sequence ambiguity.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> euler_characteristic(Z, structure_sheaf(X))
0
```
"""
function euler_characteristic(Z::ZeroLocus, F::CompletelyReducibleBundle)
  marked_dynkin_type(variety(F)) == marked_dynkin_type(Z.ambient) || throw(
    ArgumentError("euler_characteristic requires a bundle on the ambient variety.")
  )
  _euler_characteristic_from_counts(Z, _to_counts(F))
end

"""
    euler_characteristic(Z::ZeroLocus) -> BigInt

Topological Euler characteristic ``\\chi_{\\mathrm{top}}(Z) = \\sum_p (-1)^p \\chi(Z, \\Omega^p_Z)``,
the alternating sum of the Euler characteristics of the exterior powers of the
cotangent bundle. Each term is exact (no long-exact-sequence ambiguity), so the
result is determined even when the Hodge diamond is not, and it is multiplicative
over products. For ``\\chi(Z, \\mathcal{O}_Z)`` use
`euler_characteristic(Z, structure_sheaf(ambient_variety(Z)))`.
"""
function euler_characteristic(Z::ZeroLocus)
  n_factors(Z) >= 2 && return prod(euler_characteristic, factors(Z))
  d = dimension(Z)
  conormal = _conormal_data(Z, structure_sheaf(Z.ambient), d)
  sum((-1)^p * _chi_row(conormal, p) for p in 0:d)
end

"""
    euler_characteristic_tangent_bundle(Z::ZeroLocus) -> BigInt

Compute the topological Euler characteristic of the tangent bundle of ``Z``,
``\\chi(Z, \\mathrm{T}_Z) = \\chi(Z, \\mathrm{T}_X|_Z) - \\chi(Z, N_{Z/X})``,
via the tangent normal sequence ``0 \\to \\mathrm{T}_Z \\to \\mathrm{T}_X|_Z \\to N_{Z/X} \\to 0``.

Both Euler characteristics are computed exactly via Koszul (no long exact
sequence ambiguity), so the result is always a precise integer.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));  # CY3 quintic

julia> euler_characteristic_tangent_bundle(Z)
-100
```
"""
function euler_characteristic_tangent_bundle(Z::ZeroLocus)
  X = ambient_variety(Z)
  euler_characteristic(Z, tangent_bundle(X)) -
  euler_characteristic(Z, defining_bundle(Z))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Cohomology of restrictions (via Koszul + LES)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cohomology_on_restriction(Z::ZeroLocus, F::CompletelyReducibleBundle)
      -> (Cohomology{BigInt}, Bool)

Compute ``\\mathrm{H}^\\bullet(Z, \\mathcal{F}|_Z)`` by:
1. Computing ``\\mathrm{H}^\\bullet(X, \\mathcal{F} \\otimes \\wedge^i \\mathcal{E}^\\vee)`` for each Koszul term via BWB.
2. Solving the Koszul filtration via the long exact sequence.

Returns `(H*(\\mathcal{F}|_Z), determined)` where `determined` indicates whether
all cohomology groups are uniquely determined by the LES.

When the Koszul filtration leaves some groups undetermined, a Serre duality
fallback is attempted: if ``\\mathrm{H}^\\bullet(Z, \\mathcal{F}^\\vee|_Z)`` is fully determined, then
``\\mathrm{H}^k(Z, \\mathcal{F}|_Z) = \\mathrm{H}^{d-k}(Z, \\mathcal{F}^\\vee|_Z)`` by Serre duality (valid when
``\\mathrm{K}_Z \\cong \\mathcal{O}_Z``, e.g. for Calabi–Yau and hyperkähler zero loci).
"""
function cohomology_on_restriction(
  Z::ZeroLocus,
  F::CompletelyReducibleBundle,
)
  marked_dynkin_type(variety(F)) == marked_dynkin_type(Z.ambient) || throw(
    ArgumentError(
      "The bundle F must live on the ambient variety of the zero locus Z"
    ),
  )
  d_Z = dimension(Z)

  koszul_cohos = _koszul_dimensions(Z, F)

  (H, determined) = solve_koszul_filtration(koszul_cohos, d_Z)
  determined && return (H, true)

  # Serre duality fallback: H^k(Z, \\mathcal{F}|_Z) = H^{d-k}(Z, F^*|_Z), valid because
  # is_calabi_yau guarantees ω_Z ≅ O_Z via adjunction.
  if is_calabi_yau(Z)
    koszul_cohos_dual = _koszul_dimensions(Z, dual(F))
    (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)
    if det_dual
      entries = BigInt[H_dual[d_Z - k] for k in 0:d_Z]
      return (Cohomology{BigInt}(entries, d_Z), true)
    end
  end

  (H, false)
end

"""
    cohomology_on_restriction(Z::ZeroLocus) -> (Cohomology{BigInt}, Bool)

Compute ``\\mathrm{H}^\\bullet(Z, \\mathcal{O}_Z)`` via the Koszul resolution.
"""
function cohomology_on_restriction(Z::ZeroLocus)
  cohomology_on_restriction(Z, structure_sheaf(Z.ambient))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Symbolic cohomology on restrictions
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cohomology_on_restriction_symbolic(
      Z::ZeroLocus, F::CompletelyReducibleBundle,
      var_counter::Ref{Int}
    ) -> Cohomology{AffineExpr}

Symbolic version of [`cohomology_on_restriction`](@ref): entries the long
exact sequences do not determine contain fresh symbolic variables instead.
"""
function cohomology_on_restriction_symbolic(
  Z::ZeroLocus,
  F::CompletelyReducibleBundle,
  var_counter::Ref{Int},
)
  entries = _restrict_to_zero_locus_les(Z, _to_counts(F), var_counter)
  Cohomology{AffineExpr}(entries, dimension(Z))
end

"""
    cohomology_on_restriction_symbolic(
      Z::ZeroLocus, var_counter::Ref{Int}
    ) -> Cohomology{AffineExpr}

Symbolic ``\\mathrm{H}^\\bullet(Z, \\mathcal{O}_Z)`` via the Koszul resolution.
"""
function cohomology_on_restriction_symbolic(
  Z::ZeroLocus,
  var_counter::Ref{Int},
)
  cohomology_on_restriction_symbolic(Z, structure_sheaf(Z.ambient), var_counter)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Alternative Koszul restriction
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _restrict_to_zero_locus_les(Z, F, var_counter) -> Vector{AffineExpr}

Compute ``\\mathrm{H}^\\bullet(Z, \\mathcal{F}|_Z)`` using the alternative LES solver.

Instead of parametrising connecting-map ranks (``\\delta``-variables),
creates a fresh symbolic variable for each output entry and applies
the alternating-sum LES equations.

Falls back to the numeric path when fully determined.
"""
function _restrict_to_zero_locus_les(
  Z::ZeroLocus, f_counts::Dict{IrrepLevi,Int}, var_counter::Ref{Int}
)
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  _restrict_to_zero_locus_les(Z, f_counts, var_counter, wedge_counts)
end

function _restrict_to_zero_locus_les(
  Z::ZeroLocus, f_counts::Dict{IrrepLevi,Int}, var_counter::Ref{Int},
  wedge_counts::Vector{Dict{IrrepLevi,Int}},
)
  d_Z = dimension(Z)

  # Compute Koszul cohomologies once (memory-efficient path)
  koszul_cohos = _koszul_dimensions(Z, f_counts, wedge_counts)

  # Try numeric solve first
  (H_numeric, det_numeric) = solve_koszul_filtration(koszul_cohos, d_Z)
  if det_numeric
    return AffineExpr[AffineExpr(H_numeric[k]) for k in 0:d_Z]
  end

  koszul_cohos_dual = nothing
  if is_calabi_yau(Z)
    # Serre duality fallback before the symbolic path (valid: ω_Z ≅ O_Z).
    koszul_cohos_dual = _koszul_dimensions(Z, _dual_counts(f_counts), wedge_counts)
    (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)
    if det_dual
      return AffineExpr[AffineExpr(H_dual[d_Z - k]) for k in 0:d_Z]
    end
  end

  # Symbolic path: chain the Koszul LES over K_r, K_{r-1}, …, K_0, harvesting
  # the exactness inequalities of every long exact sequence along the way.
  # For ω_Z ≅ O_Z the dual bundle is chained as well, and Serre duality
  # H^k(Z, \\mathcal{F}|_Z) = H^{d-k}(Z, F^∨|_Z) is imposed entry by entry: both chains
  # are sound parametrizations, so equating them is a sound constraint
  # (unlike cross-validating two undetermined numeric guesses).
  inequalities = AffineExpr[]
  primal = long_exact_sequence_cokernel(
    _reversed_koszul_vecs(koszul_cohos), var_counter; inequalities
  )
  dual_chain = if koszul_cohos_dual === nothing
    AffineExpr[]
  else
    long_exact_sequence_cokernel(
      _reversed_koszul_vecs(koszul_cohos_dual), var_counter; inequalities
    )
  end

  # One combined system, so every substitution stays synchronized between the
  # two chains and the harvested inequalities: slots 1:chain_length hold the
  # primal entries H^k(Z, \\mathcal{F}|_Z) for k = 0, …, d_X; the dual entries (if any)
  # follow; the inequalities come last.
  chain_length = length(primal)
  system = vcat(primal, dual_chain, inequalities)

  # The structural equations — vanishing above dim Z on both chains, and the
  # entrywise Serre duality between them — are re-applied whenever interval
  # propagation pins a variable, since each substitution can make a
  # previously unusable equation solvable.
  apply_structure! = function ()
    applied = false
    for k in (d_Z + 1):(chain_length - 1)
      is_zero_expr(system[k + 1]) || (applied |= _apply_equation!(system, system[k + 1]))
      if !isempty(dual_chain)
        dual_slot = chain_length + k + 1
        is_zero_expr(system[dual_slot]) ||
          (applied |= _apply_equation!(system, system[dual_slot]))
      end
    end
    if !isempty(dual_chain)
      for k in 0:d_Z
        serre_pair = system[k + 1] - system[chain_length + d_Z - k + 1]
        is_zero_expr(serre_pair) || (applied |= _apply_equation!(system, serre_pair))
      end
    end
    applied
  end

  apply_structure!()
  while _propagate_intervals!(system)
    apply_structure!() || break
  end

  system[1:(d_Z + 1)]
end

"""Extract the Koszul cohomologies as plain vectors in reversed order K_r, …, K_0."""
function _reversed_koszul_vecs(koszul_cohos::Vector{Cohomology{BigInt}})
  d_X = koszul_cohos[1].max_degree
  Vector{BigInt}[BigInt[kc[i] for i in 0:d_X] for kc in reverse(koszul_cohos)]
end

"""
    _restrict_to_zero_locus_les(Z, F::FilteredBundle, var_counter)
      -> Vector{AffineExpr}

Compute ``\\mathrm{H}^\\bullet(Z, \\mathcal{F}|_Z)`` for a filtered bundle ``\\mathcal{F}`` on the ambient variety.

Each Koszul term ``\\mathcal{F} ⊗ ∧^k \\mathcal{E}^\\vee`` is again a filtered bundle; its cohomology
on ``X`` is the abutment of the spectral sequence of the filtration
(see `_cohomology_filtered`), and the Koszul chain is then solved
with the symbolic LES solver.  When every spectral sequence visibly
degenerates, the terms are exact and the numeric solver (with its Serre
duality fallback for ``ω_Z ≅ \\mathcal{O}_Z``) is used instead.
"""
function _restrict_to_zero_locus_les(
  Z::ZeroLocus, F::FilteredBundle, var_counter::Ref{Int}
)
  marked_dynkin_type(variety(F)) == marked_dynkin_type(Z.ambient) || throw(
    ArgumentError(
      "_restrict_to_zero_locus_les requires a bundle on the ambient variety."
    ),
  )
  d_Z = dimension(Z)
  d_X = dimension(Z.ambient)
  wedges = _koszul_wedges!(Z)
  koszul_entries = [
    _cohomology_filtered(tensor_product(F, wedge), var_counter) for wedge in wedges
  ]

  if all(is_determined, Iterators.flatten(koszul_entries))
    koszul_cohos = [
      Cohomology{BigInt}(_determined_bigints(entries), d_X) for entries in koszul_entries
    ]
    (H, determined) = solve_koszul_filtration(koszul_cohos, d_Z)
    determined && return AffineExpr[AffineExpr(H[k]) for k in 0:d_Z]

    # Serre duality fallback: sound when ω_Z ≅ O_Z, provided the spectral
    # sequences of the dual Koszul terms degenerate as well.
    if is_calabi_yau(Z)
      F_dual = dual(F)
      scratch = Ref(0)
      dual_entries = [
        _cohomology_filtered(tensor_product(F_dual, wedge), scratch) for wedge in wedges
      ]
      if all(is_determined, Iterators.flatten(dual_entries))
        koszul_cohos_dual = [
          Cohomology{BigInt}(_determined_bigints(entries), d_X) for entries in dual_entries
        ]
        (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)
        det_dual && return AffineExpr[AffineExpr(H_dual[d_Z - k]) for k in 0:d_Z]
      end
    end
  end

  # Symbolic Koszul chain in reversed order K_r, …, K_0, with the exactness
  # inequalities, then vanishing H^k(Z, \\mathcal{F}|_Z) = 0 for k > dim Z and interval
  # propagation over the combined system.
  inequalities = AffineExpr[]
  chain = long_exact_sequence_cokernel(reverse(koszul_entries), var_counter; inequalities)
  system = vcat(chain, inequalities)
  for k in (d_Z + 1):d_X
    is_zero_expr(system[k + 1]) || _apply_equation!(system, system[k + 1])
  end
  _propagate_intervals!(system)
  system[1:(d_Z + 1)]
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Calabi–Yau detection
# ═══════════════════════════════════════════════════════════════════════════════

"""
    is_calabi_yau(Z::ZeroLocus) -> Bool

Check whether the zero locus ``Z`` has trivial (anti)canonical bundle:
``\\det(\\mathcal{E}) \\cong \\omega_X^\\vee``, equivalently ``c_1(Z) = 0``.

This is a necessary condition for ``Z`` to be Calabi–Yau.  For the full
Calabi–Yau check (including cohomology vanishing), see [`is_strict_calabi_yau`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> is_calabi_yau(zero_locus(line_bundle(X, 5)))
true

julia> is_calabi_yau(zero_locus(line_bundle(X, 3)))
false
```
"""
function is_calabi_yau(Z::ZeroLocus)
  dimension(Z) == 0 && return euler_characteristic(Z) == 1
  det(Z.defining_bundle) == anticanonical_bundle(Z.ambient)
end

"""
    is_strict_calabi_yau(Z::ZeroLocus) -> Bool

Check whether the zero locus ``Z`` is a strict Calabi–Yau variety:
1. ``c_1(Z) = 0`` (equivalently ``\\det(\\mathcal{E}) \\cong \\omega_X^\\vee``,
   or trivially when ``\\dim Z = 0``)
2. ``\\mathrm{H}^0(Z, \\mathcal{O}_Z) = 1`` (``Z`` is connected)
3. ``\\mathrm{H}^i(Z, \\mathcal{O}_Z) = 0`` for ``0 < i < \\dim Z``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> is_strict_calabi_yau(zero_locus(line_bundle(X, 5)))
true

julia> is_strict_calabi_yau(zero_locus(line_bundle(X, 3)))
false
```
"""
function is_strict_calabi_yau(Z::ZeroLocus)
  is_calabi_yau(Z) || return false

  d = dimension(Z)
  (H, _) = cohomology_on_restriction(Z)
  H[0] == 1 && all(H[i] == 0 for i in 1:(d - 1))
end

"""
    is_strongly_fano(Z::ZeroLocus) -> Bool

Check whether the zero locus ``Z`` is strongly Fano: the anticanonical bundle
``\\omega_Z^\\vee`` is ample on the ambient ``\\mathrm{G}/\\mathrm{P}``.

By the adjunction formula ``\\omega_Z^\\vee = (\\omega_X^\\vee \\otimes \\det(\\mathcal{E})^\\vee)|_Z``,
this holds when every Picard-basis coordinate of ``\\omega_X^\\vee \\otimes \\det(\\mathcal{E})^\\vee``
is strictly positive.  Strong Fano implies (ordinary) Fano.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> is_strongly_fano(zero_locus(line_bundle(X, 3)))  # -K_Z = O(2), ample on ℙ⁴
true

julia> is_strongly_fano(zero_locus(line_bundle(X, 5)))  # -K_Z = O(0), not ample
false
```
"""
function is_strongly_fano(Z::ZeroLocus)
  antican_Z = tensor_product(anticanonical_bundle(Z.ambient), dual(det(Z.defining_bundle)))
  all(>(0), picard_degrees(antican_Z))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Fano index of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fano_index(Z::ZeroLocus) -> Int

The Fano index of the zero locus ``Z``, defined (when ``\\mathrm{Pic}(Z) \\cong
\\mathbb{Z}``) as the unique positive integer ``r`` such that
``\\omega_Z^\\vee = r\\,H`` where ``H`` is the restriction of the ample generator of
``\\mathrm{Pic}(X)`` to ``Z``.

Computed via the adjunction formula: ``\\mathrm{K}_Z = (\\mathrm{K}_X \\otimes \\det \\mathcal{E})|_Z``, giving

```math
r_Z = r_X - \\deg(\\det \\mathcal{E}),
```

where ``r_X = \\mathop{\\mathrm{fano\\_index}}(X)`` and ``\\deg(\\det \\mathcal{E})`` is the degree of ``\\det(\\mathcal{E})`` as a multiple of the
ample generator ``\\omega_m``.

Requires the ambient variety to have Picard rank 1 (i.e., `picard_rank(ambient_variety(Z)) == 1`).
For higher-rank ambient Picard groups, use `anticanonical_degrees` and
`det` directly.

Throws an `ArgumentError` if the ambient Picard rank exceeds 1.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 3));

julia> fano_index(Z)  # -K_Z = O(5-3) = O(2)
2

julia> Z_CY = zero_locus(line_bundle(X, 5));

julia> fano_index(Z_CY)  # Calabi–Yau: -K_Z = O(0)
0
```
"""
function fano_index(Z::ZeroLocus)
  marked = marked_nodes(Z.ambient)
  length(marked) == 1 || throw(
    ArgumentError(
      "fano_index is only defined for zero loci in Picard-rank-1 ambient varieties; " *
      "use anticanonical_degrees and det for the general case."),
  )
  m = marked[1]
  det_E = det(Z.defining_bundle)
  deg_det = p_dominant_weight(only(det_E.components)).vec[m]
  fano_index(Z.ambient) - deg_det
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Euler characteristics of cotangent powers (K-theory, always exact)
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compute ``\\chi(\\Omega^p_Z)`` using the conormal recursion.

From the conormal sequence ``0 \\to \\mathcal{E}^\\vee|_Z \\to \\Omega_X|_Z \\to \\Omega_Z \\to 0``,
the K-theory relation ``[\\wedge^p \\Omega_X|_Z] = \\sum_i [\\wedge^i \\mathcal{E}^\\vee|_Z \\otimes
\\Omega^{p-i}_Z]`` gives a recursion for ``\\chi(\\Omega^p_Z)`` in terms of
Koszul-computable Euler characteristics on ``X``.
"""
function _chi_omega_p_conormal(Z::ZeroLocus, p::Int)
  _chi_omega_tensor_counts(Z, p, _to_counts(structure_sheaf(Z.ambient)))
end

# ─── Dict-based tensor and euler operations for _chi_omega_tensor ────────────

"""
Tensor product of two multiplicity dicts, returning a new multiplicity dict.
"""
function _tensor_product_counts(
  a_counts::Dict{IrrepLevi,Int}, b_counts::Dict{IrrepLevi,Int}
)
  result = Dict{IrrepLevi,Int}()
  for (a, ma) in a_counts
    for (b, mb) in b_counts
      total = ma * mb
      for (c, mult) in _tensor_product_terms(a, b)
        result[c] = get(result, c, 0) + total * mult
      end
    end
  end
  result
end

"""
Compute ``\\chi(Z, \\mathcal{F}|_Z)`` from a multiplicity dict, without creating CRBs.

Uses the Koszul spectral sequence:
``\\chi(Z, \\mathcal{F}|_Z) = \\sum_i (-1)^i \\chi(X, \\mathcal{F} \\otimes \\wedge^i \\mathcal{E}^\\vee)``
"""
function _euler_characteristic_from_counts(
  Z::ZeroLocus, f_counts::Dict{IrrepLevi,Int}
)
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  _euler_characteristic_from_counts(Z, f_counts, wedge_counts)
end

function _euler_characteristic_from_counts(
  Z::ZeroLocus, f_counts::Dict{IrrepLevi,Int},
  wedge_counts::Vector{Dict{IrrepLevi,Int}},
)
  result = BigInt(0)

  for (wi, w_counts) in enumerate(wedge_counts)
    chi = BigInt(0)
    for (a, ma) in f_counts
      for (b, mb) in w_counts
        total = ma * mb
        for (deg, dim) in _bwb_pair(a, b)
          chi += (iseven(deg) ? 1 : -1) * total * dim
        end
      end
    end
    result += (-1)^(wi - 1) * chi
  end

  result
end

"""
Recursively compute ``\\chi(Z, \\Omega^j_Z \\otimes G|_Z)`` where ``G`` is
represented as a multiplicity dict (restricted to ``Z`` via Koszul).

Works entirely with `Dict{IrrepLevi,Int}` to avoid creating large intermediate
`CompletelyReducibleBundle` objects.  Results are memoized via `memo` to avoid
exponential recomputation of shared subtrees.
"""
function _chi_omega_tensor_counts_cached(
  Z::ZeroLocus, j::Int, g_counts::Dict{IrrepLevi,Int},
  memo::Dict{Tuple{Int,UInt},BigInt},
  omegas_counts::Vector{Dict{IrrepLevi,Int}},
  wedge_counts::Vector{Dict{IrrepLevi,Int}},
  tp_memo::Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}},
)
  # Use hash of (j, g_counts) as memo key
  g_hash = hash(g_counts)
  memo_key = (j, g_hash)
  cached = get(memo, memo_key, nothing)
  cached !== nothing && return cached

  if j == 0
    result = _euler_characteristic_from_counts(Z, g_counts, wedge_counts)
    memo[memo_key] = result
    return result
  end

  ω_counts = omegas_counts[j + 1]
  r = length(wedge_counts) - 1

  fg_counts = get!(tp_memo, (g_hash, j, false)) do
    _tensor_product_counts(ω_counts, g_counts)
  end
  result = _euler_characteristic_from_counts(Z, fg_counts, wedge_counts)

  for i in 1:min(j, r)
    w_counts = wedge_counts[i + 1]
    gw_counts = get!(tp_memo, (g_hash, i, true)) do
      _tensor_product_counts(g_counts, w_counts)
    end
    result -= _chi_omega_tensor_counts_cached(
      Z, j - i, gw_counts, memo, omegas_counts, wedge_counts, tp_memo
    )
  end

  memo[memo_key] = result
  result
end

function _chi_omega_tensor_counts(
  Z::ZeroLocus, j::Int, g_counts::Dict{IrrepLevi,Int},
  memo::Dict{Tuple{Int,UInt},BigInt}=Dict{Tuple{Int,UInt},BigInt}(),
)
  omegas_counts = Dict{IrrepLevi,Int}[
    _to_counts(_cotangent_power(Z.ambient, k)) for k in 0:j
  ]
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  tp_memo = Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}()
  _chi_omega_tensor_counts_cached(
    Z, j, g_counts, memo, omegas_counts, wedge_counts, tp_memo
  )
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Tangent bundle cohomology
# ═══════════════════════════════════════════════════════════════════════════════

"""
    tangent_cohomology(Z::ZeroLocus) -> Cohomology{AffineExpr}

Compute ``\\mathrm{H}^\\bullet(Z, \\mathrm{T}_Z)`` via the normal bundle sequence
``0 \\to \\mathrm{T}_Z \\to \\mathrm{T}_X|_Z \\to \\mathcal{E}|_Z \\to 0``.

The restrictions ``\\mathrm{T}_X|_Z`` (through the spectral sequence of the height
filtration when ``\\mathrm{T}_X`` is not completely reducible) and ``\\mathcal{E}|_Z`` are
computed by the Koszul resolution, and the long exact sequence is solved
for the kernel term with [`les_kernel`](@ref).  The exact Euler
characteristic ``\\chi(\\mathrm{T}_Z) = \\chi(\\mathrm{T}_X|_Z) - \\chi(\\mathcal{E}|_Z)`` and, for a strict
Calabi–Yau, the vanishing ``\\mathrm{h}^0(\\mathrm{T}_Z) = \\mathrm{h}^{d-1,0} = 0`` are imposed.

Entries are `AffineExpr`s: exact integers where the constraints determine
them, symbolic otherwise.  ``\\mathrm{H}^1(Z, \\mathrm{T}_Z)`` is the space of first-order
deformations of ``Z``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> Z = zero_locus(line_bundle(projective_space(4), 5));  # quintic threefold

julia> tangent_cohomology(Z)
H¹ = 101
H² = 1
```
"""
function tangent_cohomology(Z::ZeroLocus)
  d_Z = dimension(Z)

  # For a product, H^q(T_Z) is the ∧¹T row of the (Künneth) polyvector
  # parallelogram. Only that row has to be determined — higher rows may stay
  # symbolic without forcing the monolithic fallback.
  if n_factors(Z) >= 2
    tangent_row = [hochschild_cohomology(Z)[1, q] for q in 0:d_Z]
    all(is_determined, tangent_row) && return Cohomology{AffineExpr}(tangent_row, d_Z)
  end

  var_counter = Ref(0)

  HT = _restrict_to_zero_locus_les(Z, filtered_tangent_bundle(Z.ambient), var_counter)
  HE = _restrict_to_zero_locus_les(Z, _to_counts(Z.defining_bundle), var_counter)
  inequalities = AffineExpr[]
  kernel = les_kernel(HT, HE, var_counter; inequalities)
  system = vcat(kernel, inequalities)

  # χ(T_Z) = χ(T_X|_Z) - χ(E|_Z) is exact from K-theory.
  chi =
    euler_characteristic(Z, tangent_bundle(Z.ambient)) -
    euler_characteristic(Z, Z.defining_bundle)
  kernel_length = length(kernel)
  _apply_equation!(
    system, _alternating_sum(@view system[1:kernel_length]) - AffineExpr(chi)
  )
  _propagate_intervals!(system)
  entries = system[1:kernel_length]

  # For a strict Calabi–Yau, T_Z ≅ Ω^{d-1}_Z ⊗ ω_Z^{-1} ≅ Ω^{d-1}_Z, and
  # h^0(Ω^{d-1}_Z) = h^{d-1,0} = h^{d-1}(O_Z) = 0.
  if d_Z >= 2 && !is_zero_expr(entries[1]) && is_strict_calabi_yau(Z)
    _apply_equation!(entries, entries[1])
  end

  _renumber_variables!(entries)
  Cohomology{AffineExpr}(entries, d_Z)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hilbert polynomial
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hilbert_polynomial(Z::ZeroLocus[, L::CompletelyReducibleBundle])
      -> Vector{Rational{BigInt}}

Compute the Hilbert polynomial ``P(t) = \\chi(Z, \\mathcal{L}^{\\otimes t}|_Z)`` of the
zero locus with respect to the polarization ``\\mathcal{L}``, as coefficients
``[a_0, a_1, \\ldots, a_d]`` with ``P(t) = \\sum_k a_k t^k``.

The polarization must be a line bundle on the ambient variety; for the
result to be a genuine Hilbert polynomial it should be ample.  Without an
explicit polarization the minimal ample line bundle
``\\mathcal{O}(1, \\ldots, 1)`` (degree one at every marked node) is used;
for Picard rank 1 this is the ample generator of ``\\mathrm{Pic}(X)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> hp = hilbert_polynomial(Z);

julia> length(hp) >= 3
true
```
"""
function hilbert_polynomial(Z::ZeroLocus, L::CompletelyReducibleBundle)
  marked_dynkin_type(variety(L)) == marked_dynkin_type(Z.ambient) || throw(
    ArgumentError("the polarization must live on the ambient variety.")
  )
  rank(L) == 1 || throw(
    ArgumentError("the polarization must be a line bundle (rank 1).")
  )

  # A degree-d polynomial is pinned down by d+1 values; the interpolation
  # is exact over the rationals.
  d = dimension(Z)
  values = Vector{Rational{BigInt}}(undef, d + 1)
  power = structure_sheaf(Z.ambient)
  for t in 0:d
    values[t + 1] = Rational{BigInt}(euler_characteristic(Z, power))
    t < d && (power = tensor_product(power, L))
  end
  _newton_interpolation(values)
end

function hilbert_polynomial(Z::ZeroLocus)
  X = Z.ambient
  hilbert_polynomial(Z, line_bundle(X, ones(Int, length(marked_nodes(X)))))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Display
# ═══════════════════════════════════════════════════════════════════════════════

function Base.show(io::IO, Z::ZeroLocus)
  X = Z.ambient
  r = codimension(Z)
  d = dimension(Z)
  print(io, "Z(s) ⊂ $X, dim = $d, codim = $r")
end
