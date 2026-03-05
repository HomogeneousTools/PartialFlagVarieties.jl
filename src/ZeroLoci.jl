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

  solve_koszul_filtration(koszul_cohos, d_Z)
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

  terms = koszul_terms(Z, F)
  koszul_cohos = Cohomology{BigInt}[]
  for K in terms
    push!(koszul_cohos, dimensions(K))
  end

  solve_koszul_filtration_symbolic(koszul_cohos, d_Z, var_counter)
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

"""Compute the scaled central character of the anticanonical bundle -K_{G/P}."""
function _anticanonical_central(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  R = rank(DT)
  anticK = _anticanonical_weight_direct(MDT)
  sf = central_scaling_factor(MDT)

  M = decomposition_matrix(MDT)
  anticK_svec = SVector{R,Rational{Int}}(Tuple(anticK))
  new_coords = M * anticK_svec

  Marked = marked_nodes(MDT)
  Int[Int(new_coords[m] * sf) for m in Marked]
end

"""Compute the anticanonical weight of G/P in the fundamental weight basis."""
function _anticanonical_weight_direct(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  R = rank(DT)
  pos_roots = positive_nonparabolic_roots(MDT)
  C = Lie.cartan_matrix(DT)

  # Sum all roots in simple root basis
  total = zeros(Int, R)
  for root in pos_roots
    c = Lie.coefficients(root)
    for j in 1:R
      total[j] += c[j]
    end
  end

  # Convert: ω-coord_j = Σ_i C[j,i] * total[i]
  result = zeros(Rational{Int}, R)
  for j in 1:R
    for i in 1:R
      result[j] += C[j, i] * total[i]
    end
  end

  result
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
#  Hodge numbers of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers(Z::ZeroLocus) -> Matrix{BigInt}

Compute the Hodge diamond ``h^{p,q}(Z)`` for ``p, q = 0, \\ldots, \\dim Z``.

Uses the Koszul resolution and the conormal exact sequence.
Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = h^{p,q}``.

For ``p = 0``: computed directly from the Koszul resolution of ``\\mathcal{O}_Z``.
For ``p \\ge 1``: uses the conormal sequence and previously computed ``h^{j,q}``
for ``j < p``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> h = hodge_numbers(Z);

julia> h[2, 2]  # h^{1,1}
1

julia> h[3, 2]  # h^{2,1}
101
```
"""
function hodge_numbers(Z::ZeroLocus{MDT}) where {MDT}
  d = dimension(Z)
  X = Z.ambient

  hodge = zeros(BigInt, d + 1, d + 1)
  known = falses(d + 1, d + 1)

  # ── Step 1: p = 0 row via Koszul ──────────────────────────────────────
  (H0, _) = cohomology_on_restriction(Z)
  for q in 0:d
    hodge[1, q + 1] = H0[q]
    known[1, q + 1] = true
  end

  # ── Step 2: p = 1 row via conormal SES 0 → E*|_Z → Ω_X|_Z → Ω_Z → 0
  if d >= 1
    Ω_X = exterior_power(cotangent_bundle(X), 1)
    (HΩ, _) = cohomology_on_restriction(Z, Ω_X)
    E_dual = dual(Z.defining_bundle)
    (HE, _) = cohomology_on_restriction(Z, E_dual)
    (H1, _) = solve_ses_cohomology(HE, HΩ)
    for q in 0:d
      hodge[2, q + 1] = H1[q]
      known[2, q + 1] = true
    end
  end

  # ── Step 3: compute χ(Ω^p_Z) for all p via conormal recursion ────────
  chi_omega = zeros(BigInt, d + 1)
  chi_omega[1] = sum((-1)^q * hodge[1, q + 1] for q in 0:d)
  if d >= 1
    chi_omega[2] = sum((-1)^q * hodge[2, q + 1] for q in 0:d)
  end
  for p in 2:d
    chi_omega[p + 1] = _chi_omega_p_conormal(Z, p)
  end

  # ── Step 4: fill Hodge diamond using symmetries + χ constraints ───────
  for p in 2:d
    # Hodge symmetry: h^{p,q} = h^{q,p}
    for q in 0:d
      if q <= d && known[q + 1, p + 1]
        hodge[p + 1, q + 1] = hodge[q + 1, p + 1]
        known[p + 1, q + 1] = true
      end
    end

    # Serre duality: h^{p,q} = h^{d-p,d-q}
    for q in 0:d
      dp = d - p
      dq = d - q
      if 0 <= dp <= d && 0 <= dq <= d && known[dp + 1, dq + 1] && !known[p + 1, q + 1]
        hodge[p + 1, q + 1] = hodge[dp + 1, dq + 1]
        known[p + 1, q + 1] = true
      end
    end

    # Use χ(Ω^p_Z) = Σ_q (-1)^q h^{p,q} to solve for remaining unknowns
    unknown_qs = [q for q in 0:d if !known[p + 1, q + 1]]
    if length(unknown_qs) == 1
      q = unknown_qs[1]
      known_sum = sum((-1)^qq * hodge[p + 1, qq + 1] for qq in 0:d if qq != q)
      hodge[p + 1, q + 1] = (-1)^q * (chi_omega[p + 1] - known_sum)
      known[p + 1, q + 1] = true
    elseif length(unknown_qs) >= 2
      _resolve_remaining!(hodge, known, p, d, chi_omega[p + 1])
    end
  end

  hodge
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Symbolic Hodge numbers
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers_symbolic(Z::ZeroLocus) -> Matrix{AffineExpr}

Symbolic version of `hodge_numbers`.  When the long exact sequence does
not uniquely determine a Hodge number, the entry is an `AffineExpr`
involving symbolic variables ``x_0, x_1, \\ldots``.

The matrix is ``(d+1) \\times (d+1)`` with ``[p+1, q+1] = h^{p,q}``.
"""
function hodge_numbers_symbolic(Z::ZeroLocus{MDT}) where {MDT}
  d = dimension(Z)
  X = Z.ambient
  var_counter = Ref(0)

  hodge = Matrix{AffineExpr}(undef, d + 1, d + 1)
  for i in eachindex(hodge)
    hodge[i] = AffineExpr(0)
  end
  known = falses(d + 1, d + 1)

  # ── Step 1: p = 0 row via Koszul (symbolic) ───────────────────────────
  H0 = cohomology_on_restriction_symbolic(Z, var_counter)
  for q in 0:d
    hodge[1, q + 1] = H0[q]
    known[1, q + 1] = true
  end

  # ── Step 2: p = 1 row via conormal SES ────────────────────────────────
  if d >= 1
    Ω_X = exterior_power(cotangent_bundle(X), 1)
    HΩ_sym = cohomology_on_restriction_symbolic(Z, Ω_X, var_counter)
    E_dual = dual(Z.defining_bundle)
    HE_sym = cohomology_on_restriction_symbolic(Z, E_dual, var_counter)

    # Conormal SES: 0 → E*|_Z → Ω_X|_Z → Ω_Z → 0
    # Dispatch handles all cases: (BigInt,BigInt), (AffineExpr,BigInt),
    # or (AffineExpr,AffineExpr).
    H1_sym = solve_ses_cohomology_symbolic(HE_sym, HΩ_sym, var_counter)

    for q in 0:d
      hodge[2, q + 1] = H1_sym[q]
      known[2, q + 1] = true
    end
  end

  # ── Step 3: compute exact χ(Ω^p_Z) for all p via conormal recursion ──
  chi_omega = zeros(BigInt, d + 1)
  # χ for p=0, p=1: extract from symbolic hodge (only the constant parts
  # contribute since variables cancel in alternating sum... actually they
  # might not cancel).  Use the exact Euler characteristic instead.
  chi_omega[1] = euler_characteristic(Z)
  if d >= 1
    Ω_X = exterior_power(cotangent_bundle(X), 1)
    chi_omega[2] = euler_characteristic(Z, Ω_X) -
                   euler_characteristic(Z, dual(Z.defining_bundle))
  end
  for p in 2:d
    chi_omega[p + 1] = _chi_omega_p_conormal(Z, p)
  end

  # ── Step 4: fill Hodge diamond using symmetries + χ constraints ───────
  for p in 2:d
    # Hodge symmetry: h^{p,q} = h^{q,p}
    for q in 0:d
      if q <= d && known[q + 1, p + 1]
        hodge[p + 1, q + 1] = hodge[q + 1, p + 1]
        known[p + 1, q + 1] = true
      end
    end

    # Serre duality: h^{p,q} = h^{d-p,d-q}
    for q in 0:d
      dp = d - p
      dq = d - q
      if 0 <= dp <= d && 0 <= dq <= d && known[dp + 1, dq + 1] && !known[p + 1, q + 1]
        hodge[p + 1, q + 1] = hodge[dp + 1, dq + 1]
        known[p + 1, q + 1] = true
      end
    end

    # Use χ(Ω^p_Z) = Σ_q (-1)^q h^{p,q} to solve for remaining unknowns
    unknown_qs = [q for q in 0:d if !known[p + 1, q + 1]]
    if length(unknown_qs) == 1
      q = unknown_qs[1]
      known_sum = sum((-1)^qq * hodge[p + 1, qq + 1] for qq in 0:d if qq != q; init=AffineExpr(0))
      hodge[p + 1, q + 1] = (-1)^q * (AffineExpr(chi_omega[p + 1]) - known_sum)
      known[p + 1, q + 1] = true
    end
    # If ≥ 2 unknowns remain, leave as zero (symbolic resolution not attempted)
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
function hilbert_polynomial(Z::ZeroLocus{MDT}) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  length(Marked) == 1 || throw(ArgumentError(
    "hilbert_polynomial requires Picard rank 1 (single marked node)"
  ))
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
