# ═══════════════════════════════════════════════════════════════════════════════
#  ZeroLoci.jl — Zero loci of sections of equivariant bundles
#
#  A ZeroLocus{MDT} represents the zero locus Z(s) of a regular section
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

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ZeroLocus{MDT}

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
mutable struct ZeroLocus{MDT<:MarkedDynkinType}
  const ambient::PartialFlagVariety{MDT}
  const defining_bundle::CompletelyReducibleBundle{MDT}
  # Exterior powers of the dual bundle: koszul_wedges[i+1] = ∧˾i E*.
  # Populated lazily on the first call that needs them; reused thereafter.
  koszul_wedges::Union{Nothing, Vector{CompletelyReducibleBundle{MDT}}}
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
function zero_locus(E::CompletelyReducibleBundle{MDT}) where {MDT}
  X = E.variety
  r = Int(rank_bundle(E))
  d = dimension(X)
  r <= d || throw(ArgumentError(
    "Bundle rank $r exceeds ambient dimension $d."
  ))
  ZeroLocus{MDT}(X, E, nothing)
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
function _koszul_wedges!(Z::ZeroLocus{MDT}) where {MDT}
  if Z.koszul_wedges === nothing
    r = Int(rank_bundle(Z.defining_bundle))
    E_dual = dual(Z.defining_bundle)
    Z.koszul_wedges = CompletelyReducibleBundle{MDT}[exterior_power(E_dual, i) for i in 0:r]
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
function koszul_terms(Z::ZeroLocus{MDT}, F::CompletelyReducibleBundle{MDT}) where {MDT}
  CompletelyReducibleBundle{MDT}[tensor_product(F, w) for w in _koszul_wedges!(Z)]
end

"""
    koszul_terms(Z::ZeroLocus) -> Vector{CompletelyReducibleBundle}

Return the terms of the (untwisted) Koszul complex for ``\\mathcal{O}_Z``.
"""
function koszul_terms(Z::ZeroLocus{MDT}) where {MDT}
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
function euler_characteristic(Z::ZeroLocus{MDT}, F::CompletelyReducibleBundle{MDT}) where {MDT}
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
function euler_characteristic(Z::ZeroLocus{MDT}) where {MDT}
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
  Z::ZeroLocus{MDT},
  F::CompletelyReducibleBundle{MDT},
) where {MDT}

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
      ((-1)^(i - 1)) * sum((-1)^k * koszul_cohos[i][k] for k in 0:dim_ambient; init=BigInt(0))
      for i in 1:length(koszul_cohos);
      init=BigInt(0),
    )
    chi_numeric = sum((-1)^k * H[k] for k in 0:d_Z; init=BigInt(0))

    dim_ambient_dual = koszul_cohos_dual[1].dim_variety
    chi_exact_dual = sum(
      ((-1)^(i - 1)) * sum((-1)^k * koszul_cohos_dual[i][k] for k in 0:dim_ambient_dual; init=BigInt(0))
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
function cohomology_on_restriction(Z::ZeroLocus{MDT}) where {MDT}
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
  Z::ZeroLocus{MDT},
  F::CompletelyReducibleBundle{MDT},
  var_counter::Ref{Int},
) where {MDT}
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
  Z::ZeroLocus{MDT},
  var_counter::Ref{Int},
) where {MDT}
  cohomology_on_restriction_symbolic(Z, structure_sheaf(Z.ambient), var_counter)
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
function is_calabi_yau_candidate(E::CompletelyReducibleBundle{MDT}) where {MDT}
  det_c1 = _determinant_central(E)
  antican_c1 = _anticanonical_central(MDT)
  det_c1 == antican_c1
end

"""Compute c₁(E) as scaled central character coordinates at the marked nodes."""
function _determinant_central(E::CompletelyReducibleBundle{MDT}) where {MDT}
  Marked = marked_nodes(MDT)
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
function _anticanonical_central(::Type{MDT}) where {MDT<:MarkedDynkinType}
  Marked = marked_nodes(MDT)
  K = length(Marked)
  sf = central_scaling_factor(MDT)
  DT = _ambient_type(MDT)
  Cinv = Lie.cartan_matrix_inverse(DT)
  degs = anticanonical_degrees(PartialFlagVariety{MDT}())

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
`det_bundle` directly.

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
  length(marked) == 1 || throw(ArgumentError(
    "fano_index is only defined for zero loci in Picard-rank-1 ambient varieties; " *
    "use anticanonical_degrees and det_bundle for the general case.")
  )
  m = marked[1]
  det_E = det_bundle(Z.defining_bundle)
  deg_det = p_dominant_weight(only(det_E.components)).vec[m]
  fano_index(Z.ambient) - deg_det
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
function hodge_numbers_symbolic(Z::ZeroLocus{MDT}) where {MDT}
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
  syms = CompletelyReducibleBundle{MDT}[symmetric_power(E_dual, k) for k in 0:half]
  omegas = CompletelyReducibleBundle{MDT}[exterior_power(cotangent_bundle(X), k) for k in 0:half]

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

    # Immediately apply Hodge symmetry to override conormal-computed values
    # that are already determined but potentially wrong (constant-part lower bounds).
    #
    # When a corner entry h^{p,0} (or h^{p,d}) is forced to a new value `target`
    # that differs from the determined constant the Koszul computation produced,
    # the alternating sum χ = Σ (-1)^q h^{p,q} changes by ±delta.
    #
    # For the MIDDLE ROW (p == d÷2, even d): the symbolic variables in that row
    # already satisfy χ = chi_exact for any variable assignment (they all cancel
    # in the alternating sum).  The χ constraint in the second loop is therefore
    # a no-op.  We must manually compensate h^{p, d÷2} to restore χ.
    #
    # For ALL OTHER ROWS (p ≠ d÷2): the symbolic variables do NOT all cancel in
    # χ, so the χ constraint in the second loop will self-correct by eliminating
    # one free variable.  No manual compensation is needed or safe.
    mid = d ÷ 2  # middle q-index (0-based)

    # h^{p,0} = h^{0,p} (Hodge symmetry)
    if p >= 1 && is_determined(hodge[1, p + 1])
      target = hodge[1, p + 1].constant
      entry = hodge[p + 1, 1]
      if is_determined(entry) && entry.constant != target
        delta = target - entry.constant  # change in h^{p,0}
        hodge[p + 1, 1] = AffineExpr(target)
        if p == half  # middle row: χ vars cancel → compensate h^{p,mid}
          # (-1)^mid * Δmid = -delta  →  Δmid = -delta / (-1)^mid
          comp = -delta * (iseven(mid) ? 1 : -1)
          e_mid = hodge[p + 1, mid + 1]
          hodge[p + 1, mid + 1] = AffineExpr(e_mid.constant + comp, e_mid.coeffs)
        end
        # For p ≠ half: let the χ constraint in the second loop handle it.
      else
        _apply_linear_constraint!(hodge, p + 1, 1, target)
      end
    end
    # h^{p,d} = h^{0,d-p} (Serre + Hodge symmetry; skipped when dp == p)
    dp = d - p
    if dp != p && dp >= 0 && is_determined(hodge[1, dp + 1])
      target = hodge[1, dp + 1].constant
      entry = hodge[p + 1, d + 1]
      if is_determined(entry) && entry.constant != target
        # p == half is impossible here (dp == p is excluded above),
        # so no middle-row compensation is needed.
        hodge[p + 1, d + 1] = AffineExpr(target)
        # The χ constraint in the second loop will absorb the imbalance.
      else
        _apply_linear_constraint!(hodge, p + 1, d + 1, target)
      end
    end
  end

  # ── Apply Hodge/Serre symmetry constraints to eliminate variables ─────
  for p in 0:half
    # Hodge symmetry: h^{p,0} = h^{0,p}
    if p >= 1 && is_determined(hodge[1, p + 1])
      _apply_linear_constraint!(hodge, p + 1, 1, hodge[1, p + 1].constant)
    end
    # Serre + Hodge: h^{p,d} = h^{d-p,d} = h^{0,d-p} (known from p=0 row)
    dp = d - p
    if dp != p && dp >= 0 && is_determined(hodge[1, dp + 1])
      _apply_linear_constraint!(hodge, p + 1, d + 1, hodge[1, dp + 1].constant)
    end
    # χ(Ω^p_Z) constraint
    chi_p = _chi_omega_p_conormal(Z, p)
    alt_sum = sum((-1)^q * hodge[p + 1, q + 1] for q in 0:d; init=AffineExpr(0))
    _apply_equation!(hodge, alt_sum - AffineExpr(chi_p))
  end

  # ── Cross-row Hodge symmetry: h^{p,q} = h^{q,p} for all computed rows ─
  # Apply h^{p,q} = h^{q,p} for p,q = 0..half to eliminate more variables
  for p in 0:half, q in 0:half
    p == q && continue
    expr = hodge[p + 1, q + 1] - hodge[q + 1, p + 1]
    _apply_equation!(hodge, expr)
  end

  # ── Middle-row Hodge symmetry constraint: h^{p,q} = h^{p,d-q} ────────
  if d % 2 == 0
    p = half
    for q in 0:(d ÷ 2 - 1)
      expr = hodge[p + 1, q + 1] - hodge[p + 1, d - q + 1]
      _apply_equation!(hodge, expr)
    end
  end

  # ── Fill rows p > ⌊d/2⌋ via Serre duality ────────────────────────────
  for p in (half + 1):d
    for q in 0:d
      hodge[p + 1, q + 1] = hodge[d - p + 1, d - q + 1]
    end
  end

  hodge
end

"""
Compute ``\\chi(\\Omega^p_Z)`` using the conormal recursion.

From the conormal sequence ``0 \\to E^*|_Z \\to \\Omega_X|_Z \\to \\Omega_Z \\to 0``,
the K-theory relation ``[\\wedge^p \\Omega_X|_Z] = \\sum_i [\\wedge^i E^*|_Z \\otimes
\\Omega^{p-i}_Z]`` gives a recursion for ``\\chi(\\Omega^p_Z)`` in terms of
Koszul-computable Euler characteristics on ``X``.
"""
function _chi_omega_p_conormal(Z::ZeroLocus{MDT}, p::Int) where {MDT}
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
  Z::ZeroLocus{MDT}, j::Int, G::CompletelyReducibleBundle{MDT},
) where {MDT}
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
  hodge::Matrix{BigInt}, known::BitMatrix, p::Int, d::Int, χ::BigInt,
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
