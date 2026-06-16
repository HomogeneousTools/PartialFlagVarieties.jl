# ═══════════════════════════════════════════════════════════════════════════════
#  FilteredHodge.jl — Hodge numbers of zero loci via the conormal spectral sequence
#
#  For Z = Z(s) ⊂ X = G/P, the cotangent bundle Ω^1_Z is filtered via the
#  conormal exact sequence:
#
#      0 → E^∨|_Z → Ω^1_X|_Z → Ω^1_Z → 0
#
#  This induces a filtration on each Ω^p_Z with associated graded pieces:
#
#      Gr_j(Ω^p_Z) = Sym^{p-j}(E^∨) ⊗ Ω^j_X |_Z,  j = 0, …, p
#
#  where Ω^j_X = ∧^j(Ω^1_X) is ITSELF a filtered bundle on X, inheriting
#  the root-height filtration of T_X (via Lemma 2.1 of arXiv:1606.04076).
#
#  The spectral sequence E₁^{j,q} = H^q(Z, Gr_j) is computed as follows:
#
#    1. Get the filtered cotangent bundle Ω^1_X = dual(filtered_tangent_bundle(X))
#       and its j-th exterior power ∧^j(Ω^1_X) as a FilteredBundle.
#
#    2. The bundle to restrict is F_{p,j} = Sym^{p-j}(E^∨) ⊗ ∧^j(Ω^1_X) ⊗ L
#       (a FilteredBundle twisted by a CompletelyReducibleBundle).
#
#    3. H^*(Z, F_{p,j}|_Z) is computed via the Koszul resolution:
#       for each term F_{p,j} ⊗ ∧^k(E^∨) (also a FilteredBundle), cohomology
#       on X is computed via the spectral sequence of the filtration, using
#       Borel–Weil–Bott on each graded piece.
#
#    4. The Koszul LES is assembled with long_exact_sequence_cokernel;
#       the outer conormal filtration is assembled with les_cokernel chaining.
#
#  Functions:
#   - filtered_hodge_numbers(Z)           : h^{p,q}(Z) via conormal spectral sequence
#   - filtered_twisted_hodge_numbers(Z, L): h^q(Ω^p_Z ⊗ L) via spectral sequence
#   - filtered_twisted_hodge_numbers(Z, j): integer-twist version (Pic rank 1)
# ═══════════════════════════════════════════════════════════════════════════════

export filtered_hodge_numbers, filtered_twisted_hodge_numbers

# ═══════════════════════════════════════════════════════════════════════════════
#  Level-1 helper: cohomology of a FilteredBundle with external var_counter
#
#  This mirrors cohomology(F::FilteredBundle) in SpectralSequence.jl, but uses
#  an external var_counter so that all symbolic variables (from the filtered
#  spectral sequence and from the outer Koszul/conormal LES chains) share a
#  single namespace.  No _reduce_cohomology_entries! or _renumber_variables!
#  here — those are done once at the very end.
# ═══════════════════════════════════════════════════════════════════════════════

function _cohomology_filtered(F::FilteredBundle, var_counter::Ref{Int})
  S = spectral_sequence(F)
  d = dimension(variety(F))
  iso = isotypical_components(S)
  entries = AffineExpr[AffineExpr(0) for _ in 0:d]

  for (weight, S_iso) in iso
    dim_w = BigInt(degree(weight))
    D = _diagonal_sums(S_iso)
    total_degs = sort(collect(keys(D)))

    # Introduce one symbolic variable per consecutive diagonal pair admitting a
    # potential differential.  Each variable represents the *total* killed
    # cohomology dim_w * delta_w (not the differential rank delta_w), so it
    # enters entries with coefficient ±1.  This guarantees _apply_equation!
    # can always eliminate it (the divisibility check: anything % 1 == 0).
    pair_vars = Dict{Tuple{Int,Int},AffineExpr}()
    for k in 1:(length(total_degs) - 1)
      i, j = total_degs[k], total_degs[k + 1]
      j == i + 1 || continue
      _has_potential_diff(S_iso, i) || continue
      pair_vars[(i, j)] = symbolic_variable(var_counter[])
      var_counter[] += 1
    end

    for (i, d_i) in D
      0 <= i <= d || continue
      # Deterministic contribution: dim_w copies of d_i classes
      entries[i + 1] += AffineExpr(dim_w * d_i)
      # Variable contributions: subtract with coefficient 1
      haskey(pair_vars, (i - 1, i)) && (entries[i + 1] -= pair_vars[(i - 1, i)])
      haskey(pair_vars, (i, i + 1)) && (entries[i + 1] -= pair_vars[(i, i + 1)])
    end
  end

  entries
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Level-2 helper: Koszul restriction for a FilteredBundle
#
#  Computes H^*(Z, F|_Z) for a FilteredBundle F on X = Z.ambient via
#  the Koszul resolution.  For each term F ⊗ ∧^k(E^∨) (a FilteredBundle),
#  cohomology on X is computed by _cohomology_filtered (not plain BWB).
# ═══════════════════════════════════════════════════════════════════════════════

function _restrict_filtered_to_zero_locus(
  Z::ZeroLocus, F::FilteredBundle, var_counter::Ref{Int},
)
  wedges = _koszul_wedges!(Z)
  d_Z = dimension(Z)

  # H^*(X, F ⊗ ∧^k(E^∨)) for k = 0, …, rank(E) via filtered spectral sequence
  koszul_cohos = [_cohomology_filtered(tensor_product(F, w), var_counter) for w in wedges]

  # Assemble LES in reversed Koszul order: K_r, K_{r-1}, …, K_0
  result = long_exact_sequence_cokernel(reverse(koszul_cohos), var_counter)

  # Apply vanishing H^k(Z, F|_Z) = 0 for k > d_Z
  if length(result) > d_Z + 1
    mat = reshape(copy(result), 1, length(result))
    for k in (d_Z + 1):(length(result) - 1)
      expr = mat[1, k + 1]
      is_zero_expr(expr) && continue
      _apply_equation_in_vars!(mat, expr)
    end
    return AffineExpr[mat[1, k + 1] for k in 0:d_Z]
  end

  result[1:min(d_Z + 1, length(result))]
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Engine
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _filtered_hodge_symbolic(Z, L, var_counter) -> Matrix{AffineExpr}

Compute ``h^q(Z, \\Omega^p_Z \\otimes L|_Z)`` for all ``p, q`` using the
conormal spectral sequence with the filtered cotangent bundle.

For each ``p < d``, the E₁-page entry ``E_1^{j,q} = H^q(Z, \\mathrm{Gr}_j \\otimes L|_Z)``
is computed by building ``F_{p,j} = \\mathrm{Sym}^{p-j}(E^\\vee) \\otimes
\\wedge^j(\\Omega^1_X)_{\\mathrm{filt}} \\otimes L`` as a `FilteredBundle` and
restricting via `_restrict_filtered_to_zero_locus`. For ``p = d``, ``\\omega_Z \\otimes L``
is an ordinary line bundle handled by the plain Koszul LES.
"""
function _filtered_hodge_symbolic(
  Z::ZeroLocus, L::CompletelyReducibleBundle, var_counter::Ref{Int},
)
  X = Z.ambient
  E = Z.defining_bundle
  E_dual = dual(E)
  d = dimension(Z)

  # ── Filtered cotangent bundle and its exterior powers ─────────────────────
  # Ω^1_X filtered by root height; ∧^j(Ω^1_X)_filt for j = 0..d
  Omega_filt = dual(filtered_tangent_bundle(X))
  omega_filt_powers = FilteredBundle[exterior_power(Omega_filt, j) for j in 0:d]

  # ── Precomputed data for χ constraints and p = d ─────────────────────────
  # χ uses the total bundle (forgetting filtration) — same alternating sum.
  omegas_counts = Dict{IrrepLevi,Int}[_to_counts(_cotangent_power(X, j)) for j in 0:d]
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  l_counts = _to_counts(L)

  M = fill(AffineExpr(0), d + 1, d + 1)

  for p in 0:d
    if p == d
      # ω_Z ⊗ L = (K_X ⊗ det(E)) ⊗ L is an ordinary line bundle — plain Koszul
      omega_Z_counts = _to_counts(tensor_product(canonical_bundle(X), det(E)))
      f_counts = _tensor_product_counts(omega_Z_counts, l_counts)
      Hp = _restrict_to_zero_locus_les(Z, f_counts, var_counter, wedge_counts)
    else
      # ── E₁ page: H^*(Z, Gr_j ⊗ L|_Z) for j = 0, …, p ───────────────────
      # Gr_j = Sym^{p-j}(E^∨) ⊗ ∧^j(Ω^1_X)_filt, twisted by L.
      # Each Gr_j is a FilteredBundle; its restriction to Z uses the
      # FilteredBundle spectral sequence inside the Koszul LES.
      e1 = Vector{Vector{AffineExpr}}(undef, p + 1)
      for j in 0:p
        # Sym^{p-j}(E^∨) ⊗ L: ordinary CompletelyReducibleBundle
        sym_L = tensor_product(symmetric_power(E_dual, p - j), L)
        # Twist the filtered cotangent power by sym_L
        Fj = tensor_product(omega_filt_powers[j + 1], sym_L)
        # H^*(Z, Fj|_Z) via Koszul + FilteredBundle spectral sequence
        Hj = _restrict_filtered_to_zero_locus(Z, Fj, var_counter)
        e1[j + 1] = AffineExpr[k <= length(Hj) ? Hj[k] : AffineExpr(0) for k in 1:(d + 1)]
      end

      # ── Assemble conormal spectral sequence via les_cokernel chain ────────
      Hp = e1[1]
      for k in 2:(p + 1)
        Hp = les_cokernel(Hp, e1[k], var_counter)
      end
    end

    length(Hp) > d + 1 && (Hp = Hp[1:(d + 1)])
    for q in 0:d
      M[p + 1, q + 1] = q < length(Hp) ? Hp[q + 1] : AffineExpr(0)
    end
  end

  # ── Exact boundary injection p = 0 ───────────────────────────────────────
  (H_L, _) = cohomology_on_restriction(Z, L)
  for q in 0:d
    val = BigInt(q <= H_L.dim_variety ? H_L[q] : 0)
    _inject_exact!(M, 1, q + 1, val)
  end

  # ── Exact boundary injection p = d ───────────────────────────────────────
  omega_Z_L = tensor_product(tensor_product(canonical_bundle(X), det(E)), L)
  (H_top, _) = cohomology_on_restriction(Z, omega_Z_L)
  for q in 0:d
    val = BigInt(q <= H_top.dim_variety ? H_top[q] : 0)
    _inject_exact!(M, d + 1, q + 1, val)
  end

  # ── χ constraint per row ──────────────────────────────────────────────────
  chi_memo = Dict{Tuple{Int,UInt},BigInt}()
  chi_tp_memo = Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}()
  for p in 0:d
    chi_p = _chi_omega_tensor_counts_cached(
      Z, p, l_counts, chi_memo, omegas_counts, wedge_counts, chi_tp_memo,
    )
    alt_sum = _alternating_sum(M, p + 1, d)
    eq = alt_sum - AffineExpr(chi_p)
    !isempty(eq.coeffs) && _apply_equation!(M, eq)
  end

  M
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Public API
# ═══════════════════════════════════════════════════════════════════════════════

"""
    filtered_hodge_numbers(Z::ZeroLocus) -> Matrix{AffineExpr}

Compute the Hodge numbers ``h^{p,q}(Z) = \\dim H^q(Z, \\Omega^p_Z)`` using the
conormal spectral sequence of the **filtered** cotangent bundle.

For each ``p``, ``\\Omega^p_Z`` inherits a filtration from the conormal exact
sequence ``0 \\to E^\\vee|_Z \\to \\Omega^1_X|_Z \\to \\Omega^1_Z \\to 0``,
with associated graded pieces:

```math
\\mathrm{Gr}_j(\\Omega^p_Z) = \\mathrm{Sym}^{p-j}(E^\\vee) \\otimes \\Omega^j_X \\big|_Z,
\\quad j = 0, \\ldots, p
```

Here ``\\Omega^j_X = \\bigwedge^j(\\Omega^1_X)`` is itself a **filtered bundle**
from the root-height filtration of ``T_{G/P}`` (arXiv:1606.04076 Lemma 2.1).
The E₁-page cohomology ``H^q(Z, \\mathrm{Gr}_j)`` is computed by applying the
FilteredBundle spectral sequence inside the Koszul restriction to ``Z``.

After assembling the spectral sequence, the same constraint propagation as
[`hodge_numbers`](@ref) is applied: Hodge symmetry ``h^{p,q} = h^{q,p}``,
Serre duality ``h^{p,q} = h^{d-p,d-q}``, boundary and Euler characteristic
constraints. Rows ``p > \\lfloor d/2 \\rfloor`` are filled by Serre duality.

Returns a ``(d+1) \\times (d+1)`` matrix of `AffineExpr` with
``M[p+1,\\, q+1] = h^{p,q}(Z)``. Fully determined entries display as integers;
use [`is_determined`](@ref) to check individual entries.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> M = filtered_hodge_numbers(Z);

julia> M[2, 2].constant  # h^{1,1}
1

julia> M[4, 1].constant  # h^{3,0} = 1 (Calabi-Yau)
1
```
"""
function filtered_hodge_numbers(Z::ZeroLocus)
  X = Z.ambient
  L = structure_sheaf(X)
  M = _filtered_hodge_symbolic(Z, L, Ref(0))

  d = dimension(Z)
  half = d ÷ 2

  g_counts = _to_counts(L)
  omegas_counts_half = Dict{IrrepLevi,Int}[
    _to_counts(_cotangent_power(X, k)) for k in 0:half
  ]
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  chi_memo = Dict{Tuple{Int,UInt},BigInt}()
  chi_tp_memo = Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}()
  chi_vals = BigInt[
    _chi_omega_tensor_counts_cached(
      Z, p, g_counts, chi_memo, omegas_counts_half, wedge_counts, chi_tp_memo,
    ) for p in 0:half
  ]

  constraint_changed = true
  while constraint_changed
    constraint_changed = false

    for p in 1:half
      if is_determined(M[1, p + 1])
        constraint_changed =
          _apply_hodge_constraint!(M, p + 1, 1, M[1, p + 1].constant, chi_vals[p + 1], d) ||
          constraint_changed
      end
      dp = d - p
      if dp != p && dp >= 0 && is_determined(M[1, dp + 1])
        constraint_changed =
          _apply_hodge_constraint!(M, p + 1, d + 1, M[1, dp + 1].constant, chi_vals[p + 1], d) ||
          constraint_changed
      end
    end

    for p in 0:half
      alt_sum = _alternating_sum(M, p + 1, d)
      constraint_changed =
        _apply_equation!(M, alt_sum - AffineExpr(chi_vals[p + 1])) ||
        constraint_changed
    end

    for p in 0:half, q in 0:half
      p == q && continue
      constraint_changed =
        _apply_hodge_pair!(M, p + 1, q + 1, q + 1, p + 1, chi_vals, d) ||
        constraint_changed
    end

    if d % 2 == 0
      p = half
      for q in 0:(d ÷ 2 - 1)
        constraint_changed =
          _apply_hodge_pair!(M, p + 1, q + 1, p + 1, d - q + 1, chi_vals, d) ||
          constraint_changed
      end
    end

    for p in 0:half, q in 0:d
      dq = d - q
      dp = d - p
      (0 <= dq <= half) || continue
      (0 <= dp <= d) || continue
      (p == dq && q == dp) && continue
      constraint_changed =
        _apply_hodge_pair!(M, p + 1, q + 1, dq + 1, dp + 1, chi_vals, d) ||
        constraint_changed
    end
  end

  for p in (half + 1):d
    for q in 0:d
      M[p + 1, q + 1] = M[d - p + 1, d - q + 1]
    end
  end

  _renumber_variables!(M)
  M
end

"""
    filtered_twisted_hodge_numbers(Z::ZeroLocus, L::CompletelyReducibleBundle)
        -> Matrix{AffineExpr}

Compute the twisted Hodge numbers ``h^q(Z, \\Omega^p_Z \\otimes L|_Z)``
using the conormal spectral sequence with the filtered cotangent bundle.

Each graded piece ``\\mathrm{Gr}_j = \\mathrm{Sym}^{p-j}(E^\\vee) \\otimes
\\wedge^j(\\Omega^1_X)_{\\mathrm{filt}} \\otimes L`` is a `FilteredBundle`;
its restriction to ``Z`` uses `_restrict_filtered_to_zero_locus`.

Returns a ``(d+1) \\times (d+1)`` matrix of `AffineExpr` entries where
``M[p+1,\\, q+1] = h^q(Z, \\Omega^p_Z \\otimes L|_Z)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> M = filtered_twisted_hodge_numbers(Z, structure_sheaf(X));

julia> M[1, 1].constant  # h^0(Z, 𝒪_Z)
1
```
"""
function filtered_twisted_hodge_numbers(
  Z::ZeroLocus, L::CompletelyReducibleBundle,
)
  M = _filtered_hodge_symbolic(Z, L, Ref(0))
  _renumber_variables!(M)
  M
end

"""
    filtered_twisted_hodge_numbers(Z::ZeroLocus, j::Integer) -> Matrix{AffineExpr}

Compute the twisted Hodge numbers ``h^q(Z, \\Omega^p_Z(j))`` via the conormal
spectral sequence, where the twist is ``\\mathcal{O}_X(j)|_Z``.

Requires the ambient variety to have Picard rank 1.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> M = filtered_twisted_hodge_numbers(Z, 0);

julia> M[1, 1].constant  # h^0(Z, 𝒪_Z)
1
```
"""
function filtered_twisted_hodge_numbers(Z::ZeroLocus, j::Integer)
  X = Z.ambient
  length(marked_nodes(X)) == 1 || throw(ArgumentError(
    "filtered_twisted_hodge_numbers with integer twist requires Picard rank 1",
  ))
  filtered_twisted_hodge_numbers(Z, line_bundle(X, Int(j)))
end
