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
export codimension, normal_bundle, conormal_bundle
export koszul_terms, cohomology_on_restriction, cohomology_on_restriction_symbolic
export is_calabi_yau, is_calabi_yau_candidate
export is_weak_fano
export fano_index
export hilbert_polynomial
export hodge_numbers_symbolic
export hodge_numbers_les
export euler_characteristic_tangent_bundle

"""
Renumber variables in an `AffineExpr` matrix to use contiguous IDs
starting from 0.  This makes the output more readable (e.g. `x_0, x_1, …`)
rather than reflecting the internal counter."""
function _renumber_variables!(M::Matrix{AffineExpr})
  # Collect all variable IDs
  old_ids = Set{Int}()
  for e in M
    for v in keys(e.coeffs)
      push!(old_ids, v)
    end
  end
  isempty(old_ids) && return M

  # Build mapping: sorted old IDs → 0, 1, 2, …
  mapping = Dict{Int,Int}()
  for (new_id, old_id) in enumerate(sort!(collect(old_ids)))
    mapping[old_id] = new_id - 1  # 0-based
  end

  # Apply mapping
  for i in eachindex(M)
    e = M[i]
    isempty(e.coeffs) && continue
    new_coeffs = Dict{Int,BigInt}(mapping[k] => v for (k, v) in e.coeffs)
    M[i] = AffineExpr(e.constant, new_coeffs)
  end
  M
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ZeroLocus

The zero locus ``Z(s)`` of a regular section ``s \\in H^0(X, E)`` of an
equivariant bundle ``E`` on the partial flag variety ``X = G/P``.

Assumes the section is regular, so ``\\dim Z = \\dim X - \\mathrm{rank}(E)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> dimension(Z)
3

julia> euler_characteristic(Z)
0
```
"""
mutable struct ZeroLocus
  const ambient::PartialFlagVariety
  const defining_bundle::CompletelyReducibleBundle
  # Exterior powers of the dual bundle: koszul_wedges[i+1] = ∧˾i E*.
  # Populated lazily on the first call that needs them; reused thereafter.
  koszul_wedges::Union{Nothing,Vector{CompletelyReducibleBundle}}
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Constructors
# ═══════════════════════════════════════════════════════════════════════════════

"""
    zero_locus(E::CompletelyReducibleBundle) -> ZeroLocus

Construct the zero locus of a regular section of the equivariant
bundle ``E``.  Requires ``\\mathrm{rank}(E) \\le \\dim(X)``.

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
  r = Int(rank_bundle(E))
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
  E2 = _lift_bundle_to_product(X, defining_bundle(Z2), rank(X1))
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
#  Accessors
# ═══════════════════════════════════════════════════════════════════════════════

"""
Return (and lazily compute) the cached vector `[∧°E*, ∧¹E*, …, ∧ʳE*]`.
All Koszul computations go through this accessor so the exterior powers
are computed at most once per `ZeroLocus`, regardless of how many times
different twists are requested.
"""
function _koszul_wedges!(Z::ZeroLocus)
  if Z.koszul_wedges === nothing
    r = Int(rank_bundle(Z.defining_bundle))
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

Return the bundle ``E`` whose section defines the zero locus.
"""
defining_bundle(Z::ZeroLocus) = Z.defining_bundle

"""
    codimension(Z::ZeroLocus) -> Int

Return the codimension of ``Z`` in ``X``, equal to ``\\mathrm{rank}(E)``.
"""
codimension(Z::ZeroLocus)::Int = rank_bundle(Z.defining_bundle)

"""
    dimension(Z::ZeroLocus) -> Int

Return the dimension of the zero locus ``Z = \\dim X - \\mathrm{rank}(E)``.
"""
dimension(Z::ZeroLocus)::Int = dimension(Z.ambient) - codimension(Z)

"""
    normal_bundle(Z::ZeroLocus) -> CompletelyReducibleBundle

The normal bundle ``N_{Z/X} \\cong E|_Z``.
"""
normal_bundle(Z::ZeroLocus) = Z.defining_bundle

"""
    conormal_bundle(Z::ZeroLocus) -> CompletelyReducibleBundle

The conormal bundle ``N^*_{Z/X} \\cong E^*|_Z``.
"""
conormal_bundle(Z::ZeroLocus) = dual(Z.defining_bundle)

# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul complex construction
# ═══════════════════════════════════════════════════════════════════════════════

"""
    koszul_terms(Z::ZeroLocus, F::CompletelyReducibleBundle)
      -> Vector{CompletelyReducibleBundle}

Return the terms of the twisted Koszul complex:
``[F \\otimes \\wedge^0 E^*, F \\otimes \\wedge^1 E^*, \\ldots,
   F \\otimes \\wedge^r E^*]``
where ``E`` is the defining bundle and ``r = \\mathrm{rank}(E)``.
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
    maxsize=_cache_maxsize(b, _DEFAULT_BWB_FRAC),
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

"""
Compute and cache the BWB contributions `[(deg, dim), ...]` for
`tensor_product(a, b)`.  Returns an empty vector when all components are
acyclic.
"""
function _bwb_pair(a::IrrepLevi, b::IrrepLevi)
  get!(_BWB_PAIR_CACHE, (a, b)) do
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

Compute dimension-valued cohomology for each Koszul term ``F ⊗ ∧^i E^*``
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

  d = dimension(F.variety)

  # Deduplicate F components once
  f_counts = Dict{IrrepLevi,Int}()
  for c in F.components
    f_counts[c] = get(f_counts, c, 0) + 1
  end

  wedges = _koszul_wedges!(Z)
  result = Vector{Cohomology{BigInt}}(undef, length(wedges))

  for (wi, w) in enumerate(wedges)
    # Deduplicate wedge components
    w_counts = Dict{IrrepLevi,Int}()
    for c in w.components
      w_counts[c] = get(w_counts, c, 0) + 1
    end

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
  for (c, m) in f_counts
    d = dual(c)
    result[d] = get(result, d, 0) + m
  end
  result
end

"""
    _koszul_euler_characteristics(Z::ZeroLocus, F::CompletelyReducibleBundle)
      -> Vector{BigInt}

Compute ``χ(X, F ⊗ ∧^i E^*)`` for each Koszul term without materialising
intermediate `CompletelyReducibleBundle` objects.
"""
function _koszul_euler_characteristics(Z::ZeroLocus, F::CompletelyReducibleBundle)
  marked_dynkin_type(variety(F)) == marked_dynkin_type(Z.ambient) || throw(
    ArgumentError(
      "_koszul_euler_characteristics requires a bundle on the ambient variety."
    ),
  )

  # Deduplicate F components once
  f_counts = Dict{IrrepLevi,Int}()
  for c in F.components
    f_counts[c] = get(f_counts, c, 0) + 1
  end

  wedges = _koszul_wedges!(Z)
  result = Vector{BigInt}(undef, length(wedges))

  for (wi, w) in enumerate(wedges)
    w_counts = _to_counts(w)

    chi = BigInt(0)
    for (a, ma) in f_counts
      for (b, mb) in w_counts
        total = ma * mb
        for (deg, dim) in _bwb_pair(a, b)
          chi += (iseven(deg) ? 1 : -1) * total * dim
        end
      end
    end
    result[wi] = chi
  end

  result
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Euler characteristic (always exact)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    euler_characteristic(Z::ZeroLocus, F::CompletelyReducibleBundle) -> BigInt

Compute ``\\chi(Z, F|_Z) = \\sum_{i=0}^{r} (-1)^i \\chi(X, F \\otimes \\wedge^i E^*)``.
This is always exact — no long exact sequence ambiguity.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> euler_characteristic(Z)
0
```
"""
function euler_characteristic(Z::ZeroLocus, F::CompletelyReducibleBundle)
  chis = _koszul_euler_characteristics(Z, F)
  result = BigInt(0)
  for (i, χ) in enumerate(chis)
    result += (-1)^(i - 1) * χ
  end
  result
end

"""
    euler_characteristic(Z::ZeroLocus) -> BigInt

Compute ``\\chi(Z, \\mathcal{O}_Z)``.
"""
function euler_characteristic(Z::ZeroLocus)
  euler_characteristic(Z, structure_sheaf(Z.ambient))
end

"""
    euler_characteristic_tangent_bundle(Z::ZeroLocus) -> BigInt

Compute the topological Euler characteristic of the tangent bundle of ``Z``,
``\\chi(Z, T_Z) = \\chi(Z, T_X|_Z) - \\chi(Z, N_{Z/X})``,
via the tangent normal sequence ``0 \\to T_Z \\to T_X|_Z \\to N_{Z/X} \\to 0``.

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
  euler_characteristic(Z, tangent_bundle(X)) - euler_characteristic(Z, defining_bundle(Z))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Cohomology of restrictions (via Koszul + LES)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cohomology_on_restriction(Z::ZeroLocus, F::CompletelyReducibleBundle)
      -> (Cohomology{BigInt}, Bool)

Compute ``H^*(Z, F|_Z)`` by:
1. Computing ``H^*(X, F \\otimes \\wedge^i E^*)`` for each Koszul term via BWB.
2. Solving the Koszul filtration via the long exact sequence.

Returns `(H*(F|_Z), determined)` where `determined` indicates whether
all cohomology groups are uniquely determined by the LES.

When the Koszul filtration leaves some groups undetermined, a Serre duality
fallback is attempted: if ``H^*(Z, F^*|_Z)`` is fully determined, then
``H^k(Z, F|_Z) = H^{d-k}(Z, F^*|_Z)`` by Serre duality (valid when
``K_Z \\cong \\mathcal{O}_Z``, e.g. for Calabi–Yau and hyperkähler zero loci).
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

  (H, det) = solve_koszul_filtration(koszul_cohos, d_Z)
  det && return (H, true)

  # Serre duality fallback: H^k(Z, F) = H^{d-k}(Z, F*) when K_Z = O_Z.
  # Try computing H*(Z, F*) directly (without further fallback to avoid recursion).
  F_dual = dual(F)
  koszul_cohos_dual = _koszul_dimensions(Z, F_dual)
  (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)

  # Serre duality is only valid H^k(F) = H^{d-k}(F*) when K_Z ≅ O_Z.
  # Guard all Serre-duality paths by checking that the zero locus is
  # a (candidate) Calabi–Yau (i.e. det(E) ≅ ω_X^{-1}).
  if is_calabi_yau_candidate(Z.defining_bundle)
    if det_dual
      # Apply Serre duality: H^k(Z, F) = H^{d-k}(Z, F*)
      entries = Vector{BigInt}(undef, d_Z + 1)
      for k in 0:d_Z
        entries[k + 1] = H_dual[d_Z - k]
      end
      return (Cohomology{BigInt}(entries, d_Z), true)
    end

    # Both F and F* are underdetermined by the Koszul filtration alone.
    # Cross-validate: if both numeric results satisfy the Euler characteristic
    # constraint AND are Serre-dual to each other, the pair is consistent with
    # all available constraints and the result can be trusted.
    chi_exact = _alternating_euler_characteristic(koszul_cohos)
    chi_numeric = euler_characteristic(H)

    chi_exact_dual = _alternating_euler_characteristic(koszul_cohos_dual)
    chi_numeric_dual = euler_characteristic(H_dual)

    serre_consistent = all(H[k] == H_dual[d_Z - k] for k in 0:d_Z)

    if chi_numeric == chi_exact && chi_numeric_dual == chi_exact_dual && serre_consistent
      return (H, true)
    end
  end

  (H, false)
end

"""
    cohomology_on_restriction(Z::ZeroLocus) -> (Cohomology{BigInt}, Bool)

Compute ``H^*(Z, \\mathcal{O}_Z)`` via the Koszul resolution.
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

Symbolic version of `cohomology_on_restriction`.  Introduces fresh
symbolic variables for undetermined connecting-map ranks.
"""
function cohomology_on_restriction_symbolic(
  Z::ZeroLocus,
  F::CompletelyReducibleBundle,
  var_counter::Ref{Int},
)
  d_Z = dimension(Z)

  # Compute Koszul cohomologies once (memory-efficient path)
  koszul_cohos = _koszul_dimensions(Z, F)

  # Try numeric solve first
  (H_numeric, det_numeric) = solve_koszul_filtration(koszul_cohos, d_Z)
  if det_numeric
    entries = AffineExpr[AffineExpr(Int(H_numeric[k])) for k in 0:d_Z]
    return Cohomology{AffineExpr}(entries, d_Z)
  end

  # Also try Serre duality fallback (numeric on the dual)
  serre_resolved = false
  H_dual_result = nothing
  if is_calabi_yau_candidate(Z.defining_bundle)
    F_dual = dual(F)
    koszul_cohos_dual = _koszul_dimensions(Z, F_dual)
    (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)
    if det_dual
      serre_resolved = true
      H_dual_result = H_dual
    end
  end

  if serre_resolved
    # H^k(F) = H^{d-k}(F*) is fully known from Serre duality
    entries = AffineExpr[AffineExpr(Int(H_dual_result[d_Z - k])) for k in 0:d_Z]
    return Cohomology{AffineExpr}(entries, d_Z)
  end

  # Numeric path underdetermined: fall back to symbolic filtration.
  # Reuse already-computed koszul_cohos.
  d_ambient = koszul_cohos[1].dim_variety
  H_sym_full = solve_koszul_filtration_symbolic(koszul_cohos, d_ambient, var_counter)

  # Apply vanishing constraints: H^k = 0 for k = d_Z+1..d_ambient.
  n_full = d_ambient + 1
  mat = Matrix{AffineExpr}(undef, 1, n_full)
  for k in 0:d_ambient
    mat[1, k + 1] = H_sym_full[k]
  end
  for k in (d_Z + 1):d_ambient
    # H^k(Z, F) = 0 for k > dim(Z)
    expr = mat[1, k + 1]
    if !is_determined(expr) || expr.constant != 0
      _apply_equation!(mat, expr - AffineExpr(BigInt(0)))
    end
  end
  # Extract the result at degrees 0..d_Z
  H_sym = Cohomology{AffineExpr}(AffineExpr[mat[1, k + 1] for k in 0:d_Z], d_Z)

  H_sym
end

"""
    cohomology_on_restriction_symbolic(
      Z::ZeroLocus, var_counter::Ref{Int}
    ) -> Cohomology{AffineExpr}

Symbolic ``H^*(Z, \\mathcal{O}_Z)`` via the Koszul resolution.
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

Compute ``H^*(Z, F|_Z)`` using the alternative LES solver.

Instead of parametrising connecting-map ranks (``\\delta``-variables),
creates a fresh symbolic variable for each output entry and applies
the alternating-sum LES equations.

Falls back to the numeric path when fully determined.
"""
function _restrict_to_zero_locus_les(
  Z::ZeroLocus, F::CompletelyReducibleBundle, var_counter::Ref{Int}
)
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  _restrict_to_zero_locus_les(Z, _to_counts(F), var_counter, wedge_counts)
end

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
  d_ambient = koszul_cohos[1].dim_variety

  # Try numeric solve first
  (H_numeric, det_numeric) = solve_koszul_filtration(koszul_cohos, d_Z)
  if det_numeric
    return AffineExpr[AffineExpr(Int(H_numeric[k])) for k in 0:d_Z]
  end

  # Try Serre duality fallback before symbolic path
  if is_calabi_yau_candidate(Z.defining_bundle)
    f_dual_counts = _dual_counts(f_counts)
    koszul_cohos_dual = _koszul_dimensions(Z, f_dual_counts, wedge_counts)
    (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)
    if det_dual
      entries = BigInt[H_dual[d_Z - k] for k in 0:d_Z]
      return AffineExpr[AffineExpr(Int(e)) for e in entries]
    end

    # Cross-validation: check if both are consistent
    chi_exact = _alternating_euler_characteristic(koszul_cohos)
    chi_numeric = euler_characteristic(H_numeric)

    chi_exact_dual = _alternating_euler_characteristic(koszul_cohos_dual)
    chi_numeric_dual = euler_characteristic(H_dual)

    serre_consistent = all(H_numeric[k] == H_dual[d_Z - k] for k in 0:d_Z)

    if chi_numeric == chi_exact && chi_numeric_dual == chi_exact_dual && serre_consistent
      return AffineExpr[AffineExpr(Int(H_numeric[k])) for k in 0:d_Z]
    end
  end

  # Symbolic path: reuse already-computed koszul_cohos
  # Extract BigInt vectors, reversed: K_r, K_{r-1}, …, K_0
  koszul_vecs = Vector{BigInt}[
    BigInt[kc[i] for i in 0:d_ambient] for kc in reverse(koszul_cohos)
  ]

  # Apply alternative LES chain
  result_full = long_exact_sequence_cokernel(koszul_vecs, var_counter)

  # Apply vanishing: H^k(Z, F|_Z) = 0 for k > d_Z
  if d_ambient > d_Z
    mat = reshape(copy(result_full), 1, length(result_full))
    for k in (d_Z + 1):d_ambient
      expr = mat[1, k + 1]
      is_zero_expr(expr) && continue
      _apply_equation_in_vars!(mat, expr)
    end
    return AffineExpr[mat[1, k + 1] for k in 0:d_Z]
  end

  result_full[1:(d_Z + 1)]
end

function _restrict_to_zero_locus_les(
  Z::ZeroLocus, var_counter::Ref{Int}
)
  _restrict_to_zero_locus_les(Z, structure_sheaf(Z.ambient), var_counter)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Calabi–Yau detection
# ═══════════════════════════════════════════════════════════════════════════════

"""
    is_calabi_yau_candidate(E::CompletelyReducibleBundle) -> Bool

Quick check: does ``\\det(E) \\cong \\omega_X^{-1}``?  This is the necessary
condition ``c_1(Z) = 0`` for a zero locus ``Z(s)`` to be Calabi–Yau.

Does not compute cohomology — only checks the determinant condition.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> is_calabi_yau_candidate(line_bundle(X, 5))
true

julia> is_calabi_yau_candidate(line_bundle(X, 3))
false
```
"""
function is_calabi_yau_candidate(E::CompletelyReducibleBundle)
  mdt = marked_dynkin_type(variety(E))
  det_c1 = _determinant_central(E)
  antican_c1 = _anticanonical_central(mdt)
  det_c1 == antican_c1
end

"""Compute c₁(E) as scaled central character coordinates at the marked nodes."""
function _determinant_central(E::CompletelyReducibleBundle)
  Marked = marked_nodes(variety(E))
  c1 = zeros(Int, length(Marked))
  for comp in components(E)
    r = Int(fiber_dimension(comp))
    for (j, _) in enumerate(Marked)
      c1[j] += r * comp.central[j]
    end
  end
  c1
end

"""Compute c₁(E) in the Picard basis (fundamental weight coordinates at marked nodes)."""
function _determinant_picard(E::CompletelyReducibleBundle)
  Marked = marked_nodes(variety(E))
  c1 = zeros(Int, length(Marked))
  for comp in components(E)
    r = Int(fiber_dimension(comp))
    v = p_dominant_weight(comp).vec
    for (j, m) in enumerate(Marked)
      c1[j] += r * v[m]
    end
  end
  c1
end

"""Compute the scaled central character of the anticanonical bundle -K_{G/P}.

Uses `anticanonical_degrees` to obtain the anticanonical ω-coordinates
at the marked nodes, then converts to the scaled central character used
internally by `IrrepLevi` (coordinates multiplied by `central_scaling_factor`).
"""
function _anticanonical_central(mdt::MarkedDynkinType)
  Marked = marked_nodes(mdt)
  K = length(Marked)
  sf = central_scaling_factor(mdt)
  DT = _ambient_type(mdt)
  Cinv = cartan_matrix_inverse(DT)
  degs = anticanonical_degrees(PartialFlagVariety(mdt))

  # central[i] = Σ_{j=1}^K  round(Int, Cinv[Marked[i], Marked[j]] * sf) * degs[j]
  # (same formula as _apply_central_ext applied to the anticanonical weight vector)
  result = Vector{Int}(undef, K)
  for i in 1:K
    total = 0
    for j in 1:K
      total += round(Int, Cinv[Marked[i], Marked[j]] * sf) * Int(degs[j])
    end
    result[i] = total
  end
  result
end

"""
Three-valued test for whether the zero locus has ample anticanonical bundle
(Fano) as seen from the ambient ``G/P``.

- `true`    — all central-coordinate differences ``(\\omega_X^{-1} - \\det E)_j > 0``:
              ``\\omega_Z^{-1}`` is ample on the ambient, hence ample on ``Z`` (Fano).
- `false`   — some difference is ``< 0``: ``\\omega_Z^{-1}`` is not nef on the ambient,
              hence ``Z`` is not weak Fano.
- `nothing` — all differences are ``\\ge 0`` but some equal zero: ``\\omega_Z^{-1}`` is
              nef but not ample on the ambient; Fano-ness of ``Z`` itself cannot be
              determined from the ambient Picard coordinates alone.
"""
function _is_fano_zero_locus(Z::ZeroLocus)::Union{Bool,Nothing}
  E = Z.defining_bundle
  mdt = marked_dynkin_type(Z.ambient)
  det_c1 = _determinant_central(E)
  antican_c1 = _anticanonical_central(mdt)
  diffs = [antican_c1[j] - det_c1[j] for j in eachindex(det_c1)]
  any(<(0), diffs) && return false
  all(>(0), diffs) ? true : nothing
end

"""
    is_calabi_yau(Z::ZeroLocus) -> Bool

Check whether the zero locus ``Z`` is a Calabi–Yau variety:
1. ``c_1(Z) = 0`` (equivalently ``\\det(E) \\cong \\omega_X^{-1}``)
2. ``H^i(Z, \\mathcal{O}_Z) = 0`` for ``0 < i < \\dim Z``
3. ``\\dim Z \\ge 2``

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
  d = dimension(Z)
  d < 2 && return false

  is_calabi_yau_candidate(Z.defining_bundle) || return false

  (H, _) = cohomology_on_restriction(Z)
  H[0] == 1 || return false
  for i in 1:(d - 1)
    H[i] == 0 || return false
  end
  true
end

"""
    is_weak_fano(Z::ZeroLocus) -> Bool

Check whether the zero locus ``Z`` is weak Fano: the anticanonical bundle
``\\omega_Z^{-1}`` is big and nef.

The anticanonical class of ``Z`` is computed via adjunction,

```math
\\omega_Z^{-1} = \\bigl(\\omega_X^{-1} \\otimes \\det(E)^{-1}\\bigr)\\big|_Z,
```

with Picard-basis coordinates ``\\lambda_Z[j] = d_j(-K_X) - c_1(E)[j]`` for
each marked node ``j``.

**Nef** is checked by ``\\lambda_Z[j] \\geq 0`` for all ``j``.  **Big** is
checked by verifying that the leading coefficient of the Hilbert polynomial

```math
P(t) = \\chi\\bigl(Z,\\, (\\omega_Z^{-1})^{\\otimes t}\\bigr)
```

is strictly positive, which is equivalent to ``(-K_Z)^{\\dim Z} > 0``.

This is strictly weaker than ampleness of ``\\omega_Z^{-1}``.  In particular,
a zero locus whose anticanonical class is nef on the ambient ``G/P`` — meaning
all Picard coordinates are non-negative but some are zero — can still be weak
Fano when the volume ``(-K_Z)^d`` is positive.

# Examples

A degree-3 hypersurface ``Z \\subset \\mathbb{P}^4`` is Fano and therefore also
weak Fano:

```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> is_weak_fano(zero_locus(line_bundle(X, 3)))
true
```

A Calabi–Yau quintic ``Z \\subset \\mathbb{P}^4`` has ``\\omega_Z^{-1}
\\cong \\mathcal{O}_Z`` (nef but volume zero, hence not big):

```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> is_weak_fano(zero_locus(line_bundle(X, 5)))
false
```

The complete flag variety ``X = A_2/B`` has ``\\omega_X^{-1} = \\mathcal{O}(2,2)``.
The zero locus of ``\\mathcal{O}(2,0)`` has ``\\omega_Z^{-1} = \\mathcal{O}(0,2)|_Z``,
which is nef but not ample on ``X`` (its Picard coordinate at node 1 is zero).
Nevertheless ``(-K_Z)^2 = 4 > 0``, so ``Z`` is weak Fano:

```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{2}, (1, 2));

julia> Z = zero_locus(line_bundle(X, [2, 0]));

julia> is_weak_fano(Z)
true
```
"""
function is_weak_fano(Z::ZeroLocus)
  d = dimension(Z)
  d < 1 && return false

  # Short-circuit: if ω_Z^{-1} is ample on the ambient G/P, Z is Fano,
  # hence also weak Fano — no further computation needed.
  _is_fano_zero_locus(Z) === true && return true

  X = Z.ambient

  # Picard-basis coordinates of ω_Z^{-1} = (ω_X^{-1} ⊗ det(E)^{-1})|_Z
  antican_Z = anticanonical_degrees(X) .- _determinant_picard(Z.defining_bundle)

  # Nef: all Picard-basis coordinates ≥ 0
  all(>=(0), antican_Z) || return false

  # Big: leading coefficient of χ(Z, (ω_Z^{-1})^⊗t) as polynomial in t is > 0
  values = Rational{BigInt}[]
  for t in 0:(d + 2)
    push!(values, Rational{BigInt}(euler_characteristic(Z, line_bundle(X, t .* antican_Z))))
  end
  coeffs = _lagrange_interpolation(values)
  length(coeffs) >= d + 1 && coeffs[d + 1] > 0
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Fano index of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fano_index(Z::ZeroLocus) -> Int

The Fano index of the zero locus ``Z``, defined (when ``\\mathrm{Pic}(Z) \\cong
\\mathbb{Z}``) as the unique positive integer ``r`` such that
``-K_Z = r\\,H`` where ``H`` is the restriction of the ample generator of
``\\mathrm{Pic}(X)`` to ``Z``.

Computed via the adjunction formula: ``K_Z = (K_X \\otimes \\det E)|_Z``, giving

```math
r_Z = r_X - \\deg(\\det E),
```

where ``r_X = \\mathop{\\mathrm{fano\\_index}}(X)`` and ``\\deg(\\det E)`` is the degree of ``\\det(E)`` as a multiple of the
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
#  Alternative Hodge numbers
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers_les(Z::ZeroLocus) -> Matrix{AffineExpr}

Compute the symbolic Hodge diamond using the alternative LES solver.

For each row ``p``, the conormal terms
``H^*(Z, \\mathrm{Sym}^{p-j}(E^*) \\otimes \\Omega^j_X|_Z)``
are computed via `_restrict_to_zero_locus_les` (alternative inner Koszul),
then chained via `long_exact_sequence_cokernel` (alternative outer conormal).

The free parameters in the output represent genuine degrees of freedom that
cannot be resolved from the Koszul complex and symmetry constraints alone.
"""
function hodge_numbers_les(Z::ZeroLocus)
  d = dimension(Z)
  X = Z.ambient
  E_dual = dual(Z.defining_bundle)
  var_counter = Ref(0)

  half = d ÷ 2
  syms_counts = Dict{IrrepLevi,Int}[
    _to_counts(symmetric_power(E_dual, k)) for k in 0:half
  ]
  omegas_counts = Dict{IrrepLevi,Int}[
    _to_counts(_cotangent_power(X, k)) for k in 0:half
  ]

  hodge = Matrix{AffineExpr}(undef, d + 1, d + 1)
  for i in eachindex(hodge)
    hodge[i] = AffineExpr(0)
  end

  # ── Compute rows p = 0..⌊d/2⌋ via alternative LES ──────────────────────
  restriction_wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  for p in 0:half
    if p == 0
      Hp = _restrict_to_zero_locus_les(
        Z, _to_counts(structure_sheaf(Z.ambient)), var_counter, restriction_wedge_counts
      )
    else
      # Conormal terms: H*(Z, Sym^{p-j}(E*) ⊗ Ω^j_X |_Z) for j = 0..p
      conormal_cohos = Vector{AffineExpr}[]
      for j in 0:p
        f_counts = _tensor_product_counts(syms_counts[p - j + 1], omegas_counts[j + 1])
        Hj = _restrict_to_zero_locus_les(Z, f_counts, var_counter, restriction_wedge_counts)
        push!(conormal_cohos, Hj)
      end
      # Outer conormal filtration via alternative LES chain
      Hp = long_exact_sequence_cokernel(conormal_cohos, var_counter)
    end
    for q in 0:d
      hodge[p + 1, q + 1] = Hp[q + 1]
    end
  end

  # ── Apply symmetry constraints (same loop as hodge_numbers_symbolic) ─
  chi_memo = Dict{Tuple{Int,UInt},BigInt}()
  g_counts = _to_counts(structure_sheaf(Z.ambient))
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  chi_tp_memo = Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}()
  chi_vals = BigInt[
    _chi_omega_p_conormal(
      Z, p, g_counts, chi_memo, omegas_counts, wedge_counts, chi_tp_memo
    ) for p in 0:half
  ]

  constraint_changed = true
  while constraint_changed
    constraint_changed = false

    # Hodge corner constraints: h^{p,0} = h^{0,p}, h^{p,d} = h^{0,d-p}
    for p in 1:half
      if is_determined(hodge[1, p + 1])
        constraint_changed =
          _apply_hodge_constraint!(
            hodge, p + 1, 1, hodge[1, p + 1].constant, chi_vals[p + 1], d
          ) ||
          constraint_changed
      end
      dp = d - p
      if dp != p && dp >= 0 && is_determined(hodge[1, dp + 1])
        constraint_changed =
          _apply_hodge_constraint!(
            hodge, p + 1, d + 1, hodge[1, dp + 1].constant, chi_vals[p + 1], d
          ) ||
          constraint_changed
      end
    end

    # χ(Ω^p_Z) constraint
    for p in 0:half
      alt_sum = _alternating_sum(hodge, p + 1, d)
      constraint_changed =
        _apply_equation!(hodge, alt_sum - AffineExpr(chi_vals[p + 1])) ||
        constraint_changed
    end

    # Cross-row Hodge symmetry: h^{p,q} = h^{q,p} for p,q ∈ 0..half
    for p in 0:half, q in 0:half
      p == q && continue
      constraint_changed =
        _apply_hodge_pair!(hodge, p + 1, q + 1, q + 1, p + 1, chi_vals, d) ||
        constraint_changed
    end

    # Middle-row Serre: h^{half,q} = h^{half,d-q}
    if d % 2 == 0
      p = half
      for q in 0:(d ÷ 2 - 1)
        constraint_changed =
          _apply_hodge_pair!(hodge, p + 1, q + 1, p + 1, d - q + 1, chi_vals, d) ||
          constraint_changed
      end
    end

    # Combined Hodge–Serre: h^{p,q} = h^{d-q,d-p}
    for p in 0:half, q in 0:d
      dq = d - q
      dp = d - p
      (0 <= dq <= half) || continue
      (0 <= dp <= d) || continue
      (p == dq && q == dp) && continue
      constraint_changed =
        _apply_hodge_pair!(hodge, p + 1, q + 1, dq + 1, dp + 1, chi_vals, d) ||
        constraint_changed
    end
  end

  # ── Fill rows p > ⌊d/2⌋ via Serre duality ────────────────────────────
  for p in (half + 1):d
    for q in 0:d
      hodge[p + 1, q + 1] = hodge[d - p + 1, d - q + 1]
    end
  end

  _renumber_variables!(hodge)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Symbolic Hodge constraint helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""
Apply `hodge[pi, qi] = target`, handling the case where the entry is
determined but contradicts the target.  When the entry is symbolic,
delegates to `_apply_equation!`.  When it's a determined wrong value,
forces the correct value and adjusts the diagonal entry of the same row
to preserve the alternating sum χ.
"""
function _apply_hodge_constraint!(
  hodge::Matrix{AffineExpr}, pi::Int, qi::Int, target::BigInt,
  chi_p::BigInt, d::Int,
)
  expr = hodge[pi, qi]
  if !is_determined(expr)
    return _apply_equation!(hodge, expr - AffineExpr(target))
  end
  expr.constant == target && return false

  # Determined but wrong: force the correct value.
  delta = target - expr.constant
  hodge[pi, qi] = AffineExpr(target)

  # Preserve the row Euler characteristic immediately, but bias the correction
  # toward the diagonal entry. This avoids injecting the shift into an arbitrary
  # symbolic off-diagonal term, which can create impossible negative Hodge
  # numbers once the remaining symmetry constraints are applied.
  q_forced = qi - 1  # 0-based column index
  p_idx = pi - 1  # 0-based
  comp_col = p_idx == q_forced ? -1 : p_idx
  if comp_col < 0
    for q in 0:d
      q == q_forced && continue
      if !is_determined(hodge[pi, q + 1])
        comp_col = q
        break
      end
    end
  end
  comp_col < 0 && (comp_col = q_forced == 0 ? d : 0)

  correction = -delta * (iseven(q_forced - comp_col) ? 1 : -1)
  e = hodge[pi, comp_col + 1]
  hodge[pi, comp_col + 1] = AffineExpr(e.constant + correction, copy(e.coeffs))
  true
end

"""
Apply `hodge[p1, q1] = hodge[p2, q2]`, handling the determined-but-wrong
case by forcing the entry with fewer variables to match the other.
"""
function _apply_hodge_pair!(
  hodge::Matrix{AffineExpr}, p1::Int, q1::Int, p2::Int, q2::Int,
  chi_vals::Vector{BigInt}, d::Int,
)
  e1 = hodge[p1, q1]
  e2 = hodge[p2, q2]
  expr = e1 - e2

  # If the equation has a variable, eliminate it normally
  if !isempty(expr.coeffs)
    return _apply_equation!(hodge, expr)
  end

  # Both sides are determined to the same value — nothing to do
  expr.constant == 0 && return false

  # Both sides are determined to different values (contradiction from
  # the Koszul solver).  Force the one with the wrong value.
  # Prefer to force the entry in the higher-indexed row (later in the
  # conormal sequence, more likely to have accumulated error).
  if p1 >= p2
    return _apply_hodge_constraint!(hodge, p1, q1, e2.constant, chi_vals[p1], d)
  else
    return _apply_hodge_constraint!(hodge, p2, q2, e1.constant, chi_vals[p2], d)
  end
end

# ═══════════════════════════════════════════════════════════════════════════════
"""
Compute ``\\chi(\\Omega^p_Z)`` using the conormal recursion.

From the conormal sequence ``0 \\to E^*|_Z \\to \\Omega_X|_Z \\to \\Omega_Z \\to 0``,
the K-theory relation ``[\\wedge^p \\Omega_X|_Z] = \\sum_i [\\wedge^i E^*|_Z \\otimes
\\Omega^{p-i}_Z]`` gives a recursion for ``\\chi(\\Omega^p_Z)`` in terms of
Koszul-computable Euler characteristics on ``X``.
"""
function _chi_omega_p_conormal(
  Z::ZeroLocus, p::Int,
  g_counts::Dict{IrrepLevi,Int}=_to_counts(structure_sheaf(Z.ambient)),
  memo::Dict{Tuple{Int,UInt},BigInt}=Dict{Tuple{Int,UInt},BigInt}(),
  omegas_counts::Vector{Dict{IrrepLevi,Int}}=Dict{IrrepLevi,Int}[
    _to_counts(_cotangent_power(Z.ambient, j)) for j in 0:p
  ],
  wedge_counts::Vector{Dict{IrrepLevi,Int}}=Dict{IrrepLevi,Int}[
    _to_counts(w) for w in _koszul_wedges!(Z)
  ],
  tp_memo::Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}=Dict{
    Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}
  }(),
)
  _chi_omega_tensor_counts_cached(
    Z, p, g_counts, memo, omegas_counts, wedge_counts, tp_memo
  )
end

# ─── Dict-based tensor and euler operations for _chi_omega_tensor ────────────

"""
Convert a `CompletelyReducibleBundle` to a `Dict{IrrepLevi,Int}` of counts.
"""
function _to_counts(E::CompletelyReducibleBundle)
  counts = Dict{IrrepLevi,Int}()
  for c in E.components
    counts[c] = get(counts, c, 0) + 1
  end
  counts
end

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
Compute ``\\chi(Z, F|_Z)`` from a multiplicity dict, without creating CRBs.

Uses the Koszul spectral sequence:
``\\chi(Z, F|_Z) = \\sum_i (-1)^i \\chi(X, F \\otimes \\wedge^i E^*)``
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

"""
Recursively compute ``\\chi(Z, \\Omega^j_Z \\otimes G|_Z)`` where ``G`` is
a bundle on ``X`` (restricted to ``Z`` via Koszul).

Base case: ``j = 0``, returns ``\\chi(Z, G|_Z)`` via Koszul.
Recursion: ``\\chi(Z, \\Omega^j_Z \\otimes G|_Z) =
            \\chi(Z, \\wedge^j \\Omega_X \\otimes G|_Z) -
            \\sum_{i=1}^{\\min(j,r)} \\chi(Z, \\Omega^{j-i}_Z \\otimes
            \\wedge^i E^* \\otimes G|_Z)``
"""
function _chi_omega_tensor(
  Z::ZeroLocus, j::Int, G::CompletelyReducibleBundle
)
  _chi_omega_tensor_counts(Z, j, _to_counts(G))
end

"""Try to resolve remaining unknowns in row `p` of the Hodge diamond."""
function _resolve_remaining!(
  hodge::Matrix{BigInt}, known::BitMatrix, p::Int, d::Int, χ::BigInt
)
  unknown_qs = [q for q in 0:d if !known[p + 1, q + 1]]

  if d - p == p
    # Middle row: Serre gives h^{p,q} = h^{p,d-q}
    # Group unknowns into pairs (q, d-q)
    handled = Set{Int}()
    for q in unknown_qs
      q in handled && continue
      dq = d - q
      if dq == q
        push!(handled, q)
      elseif dq in Set(unknown_qs) && !(dq in handled)
        # Paired unknowns: h^{p,q} = h^{p,d-q}
        # This halves the unknowns but doesn't resolve them
        push!(handled, q, dq)
      end
    end
  end

  # After pairing, if only one free parameter remains, solve from χ
  # For now, leave unsolved — the values determined by symmetry from
  # earlier rows are the most reliable.
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hilbert polynomial
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hilbert_polynomial(Z::ZeroLocus) -> Vector{Rational{BigInt}}

Compute the Hilbert polynomial ``P(t) = \\chi(Z, \\mathcal{O}_Z(t))`` as a
polynomial in ``t``.

Returns coefficients ``[a_0, a_1, \\ldots, a_d]`` so that
``P(t) = \\sum_k a_k t^k``.

For a Fano 4-fold ``Z`` of index ``i``, the anticanonical degree is
``(-K_Z)^4 = i^4 \\cdot 24 \\cdot a_4``, and ``h^0(-K_Z) = P(i)`` (by
Kodaira–Nakano vanishing).

Requires the ambient variety to have Picard rank 1.

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
function hilbert_polynomial(Z::ZeroLocus)
  X = Z.ambient
  d = dimension(Z)
  n_pts = d + 4  # degree-d polynomial needs d+1 points; extra for numerical stability
  values = Rational{BigInt}[]
  for t in 0:(n_pts - 1)
    Lt = line_bundle(X, t)
    push!(values, Rational{BigInt}(euler_characteristic(Z, Lt)))
  end
  _lagrange_interpolation(values)
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
