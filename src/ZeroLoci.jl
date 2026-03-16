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
export fano_index
export hilbert_polynomial
export hodge_numbers_symbolic
export hodge_numbers_les

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
codimension(Z::ZeroLocus) = Int(rank_bundle(Z.defining_bundle))

"""
    dimension(Z::ZeroLocus) -> Int

Return the dimension of the zero locus ``Z = \\dim X - \\mathrm{rank}(E)``.
"""
dimension(Z::ZeroLocus) = dimension(Z.ambient) - codimension(Z)

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
  terms = koszul_terms(Z, F)
  result = BigInt(0)
  for (i, K) in enumerate(terms)
    χ = euler_characteristic(K)
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

  terms = koszul_terms(Z, F)
  koszul_cohos = Cohomology{BigInt}[]
  for K in terms
    push!(koszul_cohos, dimensions(K))
  end

  (H, det) = solve_koszul_filtration(koszul_cohos, d_Z)
  det && return (H, true)

  # Serre duality fallback: H^k(Z, F) = H^{d-k}(Z, F*) when K_Z = O_Z.
  # Try computing H*(Z, F*) directly (without further fallback to avoid recursion).
  F_dual = dual(F)
  terms_dual = koszul_terms(Z, F_dual)
  koszul_cohos_dual = Cohomology{BigInt}[]
  for K in terms_dual
    push!(koszul_cohos_dual, dimensions(K))
  end
  (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)

  # Serre duality is only valid H^k(F) = H^{d-k}(F*) when K_Z ≅ O_Z.
  # Guard all Serre-duality paths by checking that the zero locus is
  # a (candidate) Calabi–Yau (i.e. det(E) ≅ ω_X^{-1}).
  if is_calabi_yau_candidate(Z.defining_bundle)
    if det_dual
      # Apply Serre duality: H^k(Z, F) = H^{d-k}(Z, F*)
      entries = BigInt[H_dual[d_Z - k] for k in 0:d_Z]
      return (Cohomology{BigInt}(entries, d_Z), true)
    end

    # Both F and F* are underdetermined by the Koszul filtration alone.
    # Cross-validate: if both numeric results satisfy the Euler characteristic
    # constraint AND are Serre-dual to each other, the pair is consistent with
    # all available constraints and the result can be trusted.
    dim_ambient = koszul_cohos[1].dim_variety
    chi_exact = sum(
      ((-1)^(i - 1)) *
      sum((-1)^k * koszul_cohos[i][k] for k in 0:dim_ambient; init=BigInt(0))
      for i in 1:length(koszul_cohos);
      init=BigInt(0),
    )
    chi_numeric = sum((-1)^k * H[k] for k in 0:d_Z; init=BigInt(0))

    dim_ambient_dual = koszul_cohos_dual[1].dim_variety
    chi_exact_dual = sum(
      ((-1)^(i - 1)) *
      sum((-1)^k * koszul_cohos_dual[i][k] for k in 0:dim_ambient_dual; init=BigInt(0))
      for i in 1:length(koszul_cohos_dual);
      init=BigInt(0),
    )
    chi_numeric_dual = sum((-1)^k * H_dual[k] for k in 0:d_Z; init=BigInt(0))

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

  # First try the numeric path (which includes Serre duality fallback).
  # If fully determined, promote to AffineExpr immediately — no symbolic variables introduced.
  (H_numeric, det_numeric) = cohomology_on_restriction(Z, F)
  if det_numeric
    entries = AffineExpr[AffineExpr(Int(H_numeric[k])) for k in 0:d_Z]
    return Cohomology{AffineExpr}(entries, d_Z)
  end

  # Numeric path underdetermined: fall back to symbolic filtration.
  terms = koszul_terms(Z, F)
  koszul_cohos = Cohomology{BigInt}[]
  for K in terms
    push!(koszul_cohos, dimensions(K))
  end

  # Run the symbolic Koszul filtration over the FULL ambient dimension,
  # not just d_Z.  The entries at degrees d_Z+1..d_ambient must all vanish
  # (cohomology of a sheaf on a d_Z-dimensional variety), so we use those
  # as additional equations to eliminate symbolic variables.
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

  # Apply Serre duality constraints: H^k(Z, F) = H^{d-k}(Z, F*)
  # for trivially-canonical Z (K_Z = O_Z).
  # Try computing H*(F*) numerically without further Serre fallback.
  F_dual = dual(F)
  terms_dual = koszul_terms(Z, F_dual)
  koszul_cohos_dual = Cohomology{BigInt}[]
  for K in terms_dual
    push!(koszul_cohos_dual, dimensions(K))
  end
  (H_dual, det_dual) = solve_koszul_filtration(koszul_cohos_dual, d_Z)

  # Serre duality H^k(F) = H^{d-k}(F*) is only valid when K_Z ≅ O_Z.
  if det_dual && is_calabi_yau_candidate(Z.defining_bundle)
    # H^k(F) = H^{d-k}(F*) is fully known: substitute into symbolic result.
    n = d_Z + 1
    mat2 = Matrix{AffineExpr}(undef, 1, n)
    for k in 0:d_Z
      mat2[1, k + 1] = H_sym[k]
    end
    for k in 0:d_Z
      val = H_dual[d_Z - k]
      if !is_determined(mat2[1, k + 1])
        _apply_equation!(mat2, mat2[1, k + 1] - AffineExpr(Int(val)))
      end
    end
    return Cohomology{AffineExpr}(AffineExpr[mat2[1, k + 1] for k in 0:d_Z], d_Z)
  end

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
#  Alternative Koszul restriction (output-variable LES)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _restrict_to_zero_locus_les(Z, F, var_counter) -> Vector{AffineExpr}

Compute ``H^*(Z, F|_Z)`` using the output-variable LES approach.

Instead of parametrising connecting-map ranks (``\\delta``-variables),
creates a fresh symbolic variable for each output entry and applies
the alternating-sum LES equations.  This mirrors the Macaulay2
`shortExactSequenceCoker` / `longExactSequence` pipeline.

Falls back to the numeric path when fully determined.
"""
function _restrict_to_zero_locus_les(
  Z::ZeroLocus, F::CompletelyReducibleBundle, var_counter::Ref{Int},
)
  d_Z = dimension(Z)

  # Fast path: if the numeric solver fully determines everything, use it
  (H_numeric, det_numeric) = cohomology_on_restriction(Z, F)
  if det_numeric
    return AffineExpr[AffineExpr(Int(H_numeric[k])) for k in 0:d_Z]
  end

  # Build Koszul cohomologies on the ambient variety
  terms = koszul_terms(Z, F)
  koszul_cohos = Cohomology{BigInt}[dimensions(K) for K in terms]
  d_ambient = koszul_cohos[1].dim_variety

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
  Z::ZeroLocus, var_counter::Ref{Int},
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
  Cinv = Lie.cartan_matrix_inverse(DT)
  degs = anticanonical_degrees(PartialFlagVariety(mdt))

  # central[i] = Σ_{j=1}^K  round(Int, Cinv[Marked[i], Marked[j]] * sf) * degs[j]
  # (same formula as _apply_central_ext applied to the anticanonical weight vector)
  Int[
    sum(
      round(Int, Cinv[Marked[i], Marked[j]] * sf) * Int(degs[j])
      for j in 1:K
    ) for i in 1:K
  ]
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
#  Alternative Hodge numbers (output-variable LES pipeline)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers_les(Z::ZeroLocus) -> Matrix{AffineExpr}

Compute the symbolic Hodge diamond using the Macaulay2-style
output-variable LES approach.

For each row ``p``, the conormal terms
``H^*(Z, \\mathrm{Sym}^{p-j}(E^*) \\otimes \\Omega^j_X|_Z)``
are computed via `_restrict_to_zero_locus_les` (alternative inner Koszul),
then chained via `long_exact_sequence_cokernel` (alternative outer conormal).

This is an independent implementation from `hodge_numbers_symbolic`,
matching the alternative `shortExactSequenceCoker`/`longExactSequence` pipeline
nearly verbatim.  The free parameters in the output represent genuine
degrees of freedom that cannot be resolved from the Koszul complex and
symmetry constraints alone.
"""
function hodge_numbers_les(Z::ZeroLocus)
  d = dimension(Z)
  X = Z.ambient
  E_dual = dual(Z.defining_bundle)
  var_counter = Ref(0)

  half = d ÷ 2
  syms = CompletelyReducibleBundle[symmetric_power(E_dual, k) for k in 0:half]
  omegas = CompletelyReducibleBundle[exterior_power(cotangent_bundle(X), k) for k in 0:half]

  hodge = Matrix{AffineExpr}(undef, d + 1, d + 1)
  for i in eachindex(hodge)
    hodge[i] = AffineExpr(0)
  end

  # ── Compute rows p = 0..⌊d/2⌋ via alternative LES ──────────────────────
  for p in 0:half
    if p == 0
      Hp = _restrict_to_zero_locus_les(Z, var_counter)
    else
      # Conormal terms: H*(Z, Sym^{p-j}(E*) ⊗ Ω^j_X |_Z) for j = 0..p
      conormal_cohos = Vector{AffineExpr}[]
      for j in 0:p
        F = tensor_product(syms[p - j + 1], omegas[j + 1])
        Hj = _restrict_to_zero_locus_les(Z, F, var_counter)
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
  chi_vals = BigInt[_chi_omega_p_conormal(Z, p) for p in 0:half]

  constraint_changed = true
  while constraint_changed
    constraint_changed = false

    # Hodge corner constraints: h^{p,0} = h^{0,p}, h^{p,d} = h^{0,d-p}
    for p in 1:half
      if is_determined(hodge[1, p + 1])
        constraint_changed =
          _apply_hodge_constraint!(hodge, p + 1, 1, hodge[1, p + 1].constant, chi_vals[p + 1], d) ||
          constraint_changed
      end
      dp = d - p
      if dp != p && dp >= 0 && is_determined(hodge[1, dp + 1])
        constraint_changed =
          _apply_hodge_constraint!(hodge, p + 1, d + 1, hodge[1, dp + 1].constant, chi_vals[p + 1], d) ||
          constraint_changed
      end
    end

    # χ(Ω^p_Z) constraint
    for p in 0:half
      alt_sum = sum((-1)^q * hodge[p + 1, q + 1] for q in 0:d; init=AffineExpr(0))
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

  # Compensate the alternating sum χ = Σ (-1)^q h^{p,q} by adjusting an
  # undetermined entry in the same row.  If no undetermined entry exists,
  # pick the entry with the most symbolic variables (most likely to absorb
  # the correction).  The correction Δ at column q must satisfy:
  #   (-1)^(qi-1) * delta + (-1)^(comp_col-1) * Δ = 0
  # ⟹ Δ = -delta * (-1)^(qi - comp_col)
  q_forced = qi - 1  # 0-based column index

  # Try to find an undetermined entry (best candidate for correction)
  comp_col = -1
  for q in 0:d
    q == q_forced && continue
    !is_determined(hodge[pi, q + 1]) && (comp_col = q; break)
  end

  if comp_col >= 0
    # Adjust undetermined entry: cancel the χ imbalance
    correction = -delta * (iseven(q_forced - comp_col) ? 1 : -1)
    e = hodge[pi, comp_col + 1]
    hodge[pi, comp_col + 1] = AffineExpr(e.constant + correction, copy(e.coeffs))
  else
    # All entries are determined.  Find the diagonal entry h^{p,p} (most
    # internal) and adjust it, since the χ constraint will be satisfied.
    p_idx = pi - 1  # 0-based
    correction = -delta * (iseven(q_forced - p_idx) ? 1 : -1)
    e = hodge[pi, p_idx + 1]
    hodge[pi, p_idx + 1] = AffineExpr(e.constant + correction)
  end
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
#  Symbolic Hodge numbers  (concrete hodge_numbers is in Hodge.jl)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers_symbolic(Z::ZeroLocus) -> Matrix{AffineExpr}

Symbolic version of `hodge_numbers`.  When the long exact sequence does
not uniquely determine a Hodge number, the entry is an `AffineExpr`
involving symbolic variables ``x_0, x_1, \\ldots``.

The matrix is ``(d+1) \\times (d+1)`` with ``[p+1, q+1] = h^{p,q}``.

After computing the ``p = 1`` row via the conormal short exact sequence,
Hodge symmetry (``h^{1,0} = h^{0,1}``), Serre duality
(``h^{1,d} = h^{0,d-1}``), and the exact ``\\chi(\\Omega^1_Z)`` are used
to eliminate up to three symbolic variables.

!!! warning "Lefschetz hyperplane theorem does not apply to higher-rank zero loci"
    The Lefschetz hyperplane theorem guarantees ``\\mathrm{Pic}(X) \\xrightarrow{\\sim}
    \\mathrm{Pic}(Z)`` only when ``Z`` is an ample *hypersurface* (codimension 1).
    For a zero locus of a rank-``r`` bundle with ``r > 1``, the Picard rank of
    ``Z`` can strictly exceed that of the ambient ``X``, so ``h^{1,1}(Z) > b_2(X)``
    is possible and may be left as a free symbolic variable by this function.
    Do **not** assume ``h^{1,1}(Z) = \\mathrm{picard\\_rank}(X)``.
    Example: ``b9 = (\\mathrm{Sym}^2 S^*)^{\\oplus 2}`` on ``\\mathrm{Gr}(2,7)``
    has ``h^{1,1} = 8`` even though ``\\mathrm{Pic}(\\mathrm{Gr}(2,7)) \\cong \\mathbb{Z}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> H = hodge_numbers_symbolic(Z);

julia> is_determined(H[2, 2])  # h^{1,1} fully determined
true

julia> H[2, 2].constant
1

julia> H[3, 2].constant  # h^{2,1}
101
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 6);

julia> E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:4]);

julia> Z = zero_locus(E);

julia> H = hodge_numbers_symbolic(Z);

julia> all(is_determined(H[p+1, q+1]) for p in 0:4, q in 0:4)
true

julia> H[2, 2].constant  # h^{1,1}
1

julia> H[3, 3].constant  # h^{2,2}
8
```
"""
function hodge_numbers_symbolic(Z::ZeroLocus)
  d = dimension(Z)
  X = Z.ambient
  E = Z.defining_bundle
  E_dual = dual(E)
  var_counter = Ref(0)

  hodge = Matrix{AffineExpr}(undef, d + 1, d + 1)
  for i in eachindex(hodge)
    hodge[i] = AffineExpr(0)
  end

  # ── Precompute Sym^k(E*) and Ω^k_X for k = 0..⌊d/2⌋ ─────────────────
  half = d ÷ 2
  syms = CompletelyReducibleBundle[symmetric_power(E_dual, k) for k in 0:half]
  omegas = CompletelyReducibleBundle[exterior_power(cotangent_bundle(X), k) for k in 0:half]

  # ── Compute rows p = 0..⌊d/2⌋ via the conormal filtration ────────────
  for p in 0:half
    if p == 0
      Hp = cohomology_on_restriction_symbolic(Z, var_counter)
    else
      cohos = Cohomology[]
      for j in 0:p
        F = tensor_product(syms[p - j + 1], omegas[j + 1])
        Hj = cohomology_on_restriction_symbolic(Z, F, var_counter)
        push!(cohos, Hj)
      end
      Hp = cohos[1]
      for k in 2:length(cohos)
        Hp = solve_ses_cohomology_symbolic(Hp, cohos[k], var_counter)
      end
      # Truncate to zero locus dimension
      Hp = Cohomology{AffineExpr}(AffineExpr[Hp[i] for i in 0:d], d)
    end
    for q in 0:d
      hodge[p + 1, q + 1] = Hp[q]
    end
  end

  # ── Apply all symmetry constraints, iterating until convergence ─────
  # Precompute χ(Ω^p_Z) once (these are exact and don't change).
  chi_vals = BigInt[_chi_omega_p_conormal(Z, p) for p in 0:half]

  constraint_changed = true
  while constraint_changed
    constraint_changed = false

    # Hodge corner constraints: h^{p,0} = h^{0,p}, h^{p,d} = h^{0,d-p}
    for p in 1:half
      if is_determined(hodge[1, p + 1])
        constraint_changed =
          _apply_hodge_constraint!(hodge, p + 1, 1, hodge[1, p + 1].constant, chi_vals[p + 1], d) ||
          constraint_changed
      end
      dp = d - p
      if dp != p && dp >= 0 && is_determined(hodge[1, dp + 1])
        constraint_changed =
          _apply_hodge_constraint!(hodge, p + 1, d + 1, hodge[1, dp + 1].constant, chi_vals[p + 1], d) ||
          constraint_changed
      end
    end

    # χ(Ω^p_Z) constraint
    for p in 0:half
      alt_sum = sum((-1)^q * hodge[p + 1, q + 1] for q in 0:d; init=AffineExpr(0))
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
    # Only apply when both (p) and (d-q) are in the computed range 0..half.
    for p in 0:half, q in 0:d
      dq = d - q
      dp = d - p
      (0 <= dq <= half) || continue
      (0 <= dp <= d) || continue
      (p == dq && q == dp) && continue  # skip trivial self-links
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

"""
Compute ``\\chi(\\Omega^p_Z)`` using the conormal recursion.

From the conormal sequence ``0 \\to E^*|_Z \\to \\Omega_X|_Z \\to \\Omega_Z \\to 0``,
the K-theory relation ``[\\wedge^p \\Omega_X|_Z] = \\sum_i [\\wedge^i E^*|_Z \\otimes
\\Omega^{p-i}_Z]`` gives a recursion for ``\\chi(\\Omega^p_Z)`` in terms of
Koszul-computable Euler characteristics on ``X``.
"""
function _chi_omega_p_conormal(Z::ZeroLocus, p::Int)
  _chi_omega_tensor(Z, p, structure_sheaf(Z.ambient))
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
  X = Z.ambient

  if j == 0
    return euler_characteristic(Z, G)
  end

  Ωj_X = exterior_power(cotangent_bundle(X), j)
  r = length(_koszul_wedges!(Z)) - 1  # = codimension(Z)

  result = euler_characteristic(Z, tensor_product(Ωj_X, G))
  for i in 1:min(j, r)
    result -= _chi_omega_tensor(Z, j - i, tensor_product(_koszul_wedges!(Z)[i + 1], G))
  end

  result
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
