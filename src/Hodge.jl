# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge.jl — Hodge numbers, twisted Hodge numbers, Hochschild cohomology,
#             and Hodge diamond printing
#
#  For partial flag varieties G/P:
#   - hodge_numbers(X): the diagonal Hodge diamond h^{p,q}(X)
#   - twisted_hodge_numbers(X, j): twisted Hodge numbers h^q(X, Ω^p(j))
#   - hochschild_cohomology(X): HKR decomposition HH^n = ⊕_{p+q=n} H^q(∧^p T)
#
#  For zero loci Z(s) ⊂ G/P:
#   - hodge_numbers(Z): AffineExpr Hodge diamond (symbolic for undetermined entries)
#   - twisted_hodge_numbers(Z, L): h^q(Z, Ω^p_Z ⊗ L|_Z)
#   - hochschild_cohomology(Z): HKR decomposition via twisted Hodge numbers
#   - print_hodge_diamond: ASCII Hodge diamond via PrettyTables
# ═══════════════════════════════════════════════════════════════════════════════

export hodge_numbers, twisted_hodge_numbers
export hochschild_cohomology, hochschild_dimension
export PolyvectorParallelogram
export print_hodge_diamond

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge numbers  h^{p,q}(X)  for partial flag varieties
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers(X::PartialFlagVariety) -> Matrix{BigInt}

Compute the Hodge numbers ``h^{p,q}(X) = \\dim H^q(X, \\Omega^p_X)``
of the partial flag variety ``X = G/P``.

Since ``G/P`` is rational (and simply connected), the Hodge diamond is
diagonal: ``h^{p,q} = 0`` for ``p \\neq q``, and ``h^{p,p} = b_{2p}``
where ``b_i`` are the Betti numbers.

Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = h^{p,q}``,
with ``d = \\dim X``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> H = hodge_numbers(X);

julia> H[1, 1]  # h^{0,0}
1

julia> H[2, 2]  # h^{1,1}
1

julia> H[1, 2]  # h^{0,1}
0
```
"""
function hodge_numbers(X::PartialFlagVariety)
  d = dimension(X)
  betti = betti_numbers(X)

  # For G/P: h^{p,q} = δ_{p,q} * b_{2p}
  # betti_numbers returns only even Betti numbers: betti[i] = b_{2(i-1)}
  H = zeros(BigInt, d + 1, d + 1)
  for p in 0:d
    if p + 1 <= length(betti)
      H[p + 1, p + 1] = betti[p + 1]  # betti[p+1] = b_{2p}
    end
  end
  H
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Twisted Hodge numbers  h^q(X, Ω^p(j))  for partial flag varieties
# ═══════════════════════════════════════════════════════════════════════════════

"""
    twisted_hodge_numbers(X::PartialFlagVariety, j::Integer) -> Matrix{BigInt}

Compute the twisted Hodge numbers ``h^q(X, \\Omega^p_X(j))`` for the
twist ``j``.

Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = h^q(X, \\Omega^p(j))``.

Requires ``X`` to be a generalized Grassmannian (Picard rank 1).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> H = twisted_hodge_numbers(X, 0);

julia> H[1, 1]  # h^0(ℙ³, 𝒪)
1
```
"""
function twisted_hodge_numbers(X::PartialFlagVariety, j::Integer)
  length(marked_nodes(X)) == 1 || throw(ArgumentError(
    "twisted_hodge_numbers requires Picard rank 1"
  ))

  d = dimension(X)
  H = zeros(BigInt, d + 1, d + 1)

  for p in 0:d
    # Ω^p(j) = ∧^p Ω ⊗ O(j)
    Omega_p = exterior_power(cotangent_bundle(X), p)
    Omega_p_j = twist(Omega_p, 1, Int(j))
    coh = dimensions(Omega_p_j)
    for q in 0:d
      H[p + 1, q + 1] = coh[q]
    end
  end

  H
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Polyvector Parallelogram — display type for Hochschild cohomology
# ═══════════════════════════════════════════════════════════════════════════════

"""
    PolyvectorParallelogram

The Hochschild–Kostant–Rosenberg decomposition of Hochschild cohomology:

``\\mathrm{HH}^n(X) = \\bigoplus_{p+q=n} H^q(X, \\bigwedge^p T_X)``

Stored as a matrix `data[p+1, q+1] = h^q(X, ∧^p T_X)`.

# Fields
- `data::Matrix{T}`: the HKR decomposition matrix (`BigInt` for varieties, `AffineExpr` for zero loci)
- `dim::Int`: the dimension of the variety
"""
struct PolyvectorParallelogram{T}
  data::Matrix{T}  # data[p+1, q+1] = h^q(X, ∧^p T_X)
  dim::Int
end

_pp_zero(::Type{BigInt}) = BigInt(0)
_pp_zero(::Type{AffineExpr}) = AffineExpr(0)

"""
    getindex(P::PolyvectorParallelogram, p::Int, q::Int)

Return ``h^q(X, \\bigwedge^p T_X)``. Uses 0-based indexing.
"""
function Base.getindex(P::PolyvectorParallelogram{T}, p::Int, q::Int) where T
  0 <= p <= P.dim || return _pp_zero(T)
  0 <= q <= P.dim || return _pp_zero(T)
  P.data[p + 1, q + 1]
end

"""
    euler_characteristic(P::PolyvectorParallelogram)

Compute the Euler characteristic of Hochschild cohomology:
``\\chi(\\mathrm{HH}^*) = \\sum_{p,q} (-1)^{p+q} h^q(X, \\bigwedge^p T_X)``
"""
function euler_characteristic(P::PolyvectorParallelogram{T}) where T
  result = _pp_zero(T)
  for p in 0:P.dim
    for q in 0:P.dim
      result += (-1)^(p + q) * P[p, q]
    end
  end
  result
end

"""
    hochschild_dimension(P::PolyvectorParallelogram, n::Int)

Compute ``\\dim \\mathrm{HH}^n(X) = \\sum_{p+q=n} h^q(X, \\bigwedge^p T_X)``.
"""
function hochschild_dimension(P::PolyvectorParallelogram{T}, n::Int) where T
  result = _pp_zero(T)
  for p in 0:P.dim
    q = n - p
    0 <= q <= P.dim || continue
    result += P[p, q]
  end
  result
end

# ─── Display as parallelogram ────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", P::PolyvectorParallelogram)
  d = P.dim

  # Width of widest entry for alignment
  w = maximum(length(string(P.data[i])) for i in eachindex(P.data))
  w = max(w, 1)

  println(io, "Polyvector parallelogram (dim = $d):")
  println(io)

  # Row n = p + q goes from 0 to 2d.
  # In each row: q ranges over max(0, n-d)..min(n, d),  p = n - q.
  # Columns are indexed by q (= ∧^q T direction); rows 0..d start
  # flush left, rows d+1..2d shift right — forming a parallelogram.
  for n in 0:(2 * d)
    indent = max(0, n - d) * (w + 2)
    print(io, " " ^ indent)

    entries = String[]
    for q in max(0, n - d):min(n, d)
      push!(entries, lpad(string(P[q, n - q]), w))
    end

    println(io, join(entries, "  "))
  end
end

function Base.show(io::IO, P::PolyvectorParallelogram)
  d = P.dim
  print(io, "PolyvectorParallelogram(dim=$d)")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hochschild cohomology of partial flag varieties via HKR
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hochschild_cohomology(X::PartialFlagVariety) -> PolyvectorParallelogram

Compute the Hochschild cohomology of ``X`` via the HKR decomposition:

``\\mathrm{HH}^n(X) = \\bigoplus_{p+q=n} H^q(X, \\bigwedge^p T_X)``

Returns a [`PolyvectorParallelogram`](@ref) encoding the full decomposition.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(2);

julia> P = hochschild_cohomology(X);

julia> P[0, 0]  # h⁰(𝒪) = 1
1

julia> P[2, 0]  # h⁰(∧²T) = h⁰(𝒪(3)) for ℙ²
10
```
"""
function hochschild_cohomology(X::PartialFlagVariety)
  d = dimension(X)
  data = zeros(BigInt, d + 1, d + 1)

  T = tangent_bundle(X)

  for p in 0:d
    wedge_p = exterior_power(T, p)
    coh = dimensions(wedge_p)
    for q in 0:d
      data[p + 1, q + 1] = coh[q]
    end
  end

  PolyvectorParallelogram(data, d)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge numbers of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers(Z::ZeroLocus) -> Matrix{AffineExpr}

Compute the Hodge diamond ``h^{p,q}(Z)`` for ``p, q = 0, \\ldots, \\dim Z``.

Uses the Koszul resolution and long exact sequences.
Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = h^{p,q}``.
Entries that are fully determined are plain integers (wrapped in `AffineExpr`);
entries that cannot be resolved from Koszul + symmetry constraints contain
symbolic variables.

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
function hodge_numbers(Z::ZeroLocus)
  hodge_numbers_les(Z)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Twisted Hodge numbers of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

"""
    twisted_hodge_numbers(Z::ZeroLocus, L::CompletelyReducibleBundle) -> Matrix{BigInt}

Compute the twisted Hodge numbers ``h^q(Z, \\Omega^p_Z \\otimes L|_Z)``
for the zero locus ``Z = Z(s) \\subset X = G/P`` and a line bundle ``L``
on ``X``.

Uses the conormal sequence to compute ``\\Omega^p_Z`` from the
ambient cotangent bundle and the conormal bundle ``E^*``:

```math
[\\wedge^p \\Omega_Z] = \\sum_{j=0}^{p} (-1)^{p-j}
  [\\mathrm{Sym}^{p-j}(E^*) \\otimes \\wedge^j \\Omega_X]
```

Each row ``p`` is computed by solving the resulting Koszul filtrations
on the ambient variety, then chaining via `solve_ses_cohomology`.

Returns a ``(d+1) \\times (d+1)`` matrix where ``d = \\dim Z`` and
entry ``[p+1, q+1] = h^q(Z, \\Omega^p_Z \\otimes L|_Z)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> M = twisted_hodge_numbers(Z, structure_sheaf(X));

julia> M[1, 1]  # h^0(Z, 𝒪_Z) = 1
1

julia> M[4, 1]  # h^0(Z, Ω³_Z) = h^{3,0} = 1 (CY3)
1
```
"""
function twisted_hodge_numbers(
  Z::ZeroLocus, L::CompletelyReducibleBundle
)
  X = Z.ambient
  E = Z.defining_bundle
  E_dual = dual(E)
  d = dimension(Z)

  # Precompute Sym^k(E*) for k = 0..d and Ω^j(X) for j = 0..d
  syms = CompletelyReducibleBundle[symmetric_power(E_dual, k) for k in 0:d]
  omegas = CompletelyReducibleBundle[exterior_power(cotangent_bundle(X), j) for j in 0:d]

  M = zeros(BigInt, d + 1, d + 1)

  for p in 0:d
    if p == d
      # Ω^d_Z = ω_Z = (ω_X ⊗ det E)|_Z, so Ω^d_Z ⊗ L = (ω_X ⊗ det(E) ⊗ L)|_Z.
      # Compute directly from Koszul, bypassing the conormal chain.
      F_direct = tensor_product(
        tensor_product(canonical_bundle(X), det(E)), L
      )
      (Hp, _) = cohomology_on_restriction(Z, F_direct)
      for q in 0:d
        M[p + 1, q + 1] = q <= Hp.dim_variety ? Hp[q] : BigInt(0)
      end
    else
      # Use numeric conormal chain solver
      cohos = Tuple{Cohomology{BigInt},Bool}[]
      for j in 0:p
        F = tensor_product(tensor_product(syms[p - j + 1], omegas[j + 1]), L)
        push!(cohos, cohomology_on_restriction(Z, F))
      end
      Hp = cohos[1][1]
      for k in 2:length(cohos)
        (Hp, _) = solve_ses_cohomology(Hp, cohos[k][1])
      end
      for q in 0:d
        M[p + 1, q + 1] = q <= Hp.dim_variety ? Hp[q] : BigInt(0)
      end
    end
  end

  M
end

"""
    twisted_hodge_numbers(Z::ZeroLocus, j::Integer) -> Matrix{BigInt}

Compute the twisted Hodge numbers ``h^q(Z, \\Omega^p_Z(j))`` where the
twist is by ``\\mathcal{O}_X(j)|_Z``.

Requires the ambient variety to have Picard rank 1.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> M = twisted_hodge_numbers(Z, 0);

julia> M[1, 1]  # h^0(Z, 𝒪_Z) = 1
1

julia> M[2, 2]  # h^1(Z, Ω¹_Z) = h^{1,1} = 1
1
```
"""
function twisted_hodge_numbers(Z::ZeroLocus, j::Integer)
  X = Z.ambient
  length(marked_nodes(X)) == 1 || throw(ArgumentError(
    "twisted_hodge_numbers with integer twist requires Picard rank 1"
  ))
  twisted_hodge_numbers(Z, line_bundle(X, Int(j)))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hochschild cohomology of zero loci via HKR
# ═══════════════════════════════════════════════════════════════════════════════

"""Check whether the zero locus has ample anticanonical bundle (Fano)."""
function _is_fano_zero_locus(Z::ZeroLocus)
  E = Z.defining_bundle
  mdt = marked_dynkin_type(Z.ambient)
  det_c1 = _determinant_central(E)
  antican_c1 = _anticanonical_central(mdt)
  all(antican_c1[j] > det_c1[j] for j in eachindex(det_c1))
end

"""
    hochschild_cohomology(Z::ZeroLocus) -> PolyvectorParallelogram{AffineExpr}

Compute the Hochschild cohomology of ``Z`` via the HKR decomposition:

``\\mathrm{HH}^n(Z) = \\bigoplus_{p+q=n} H^q(Z, \\bigwedge^p T_Z)``

Uses the identity ``H^q(Z, \\bigwedge^p T_Z) = H^q(Z, \\Omega^{d-p}_Z
\\otimes \\omega_Z^{-1})`` where ``d = \\dim Z``, reducing the computation
to twisted Hodge numbers with the anticanonical twist
``L = \\omega_X^{-1} \\otimes \\det(E)^{-1}`` lifted to the ambient variety.

When ``Z`` is Fano (``\\omega_Z^{-1}`` ample), Akizuki–Nakano vanishing
forces ``h^q(\\bigwedge^p T_Z) = 0`` for ``q > p``.

Returns a [`PolyvectorParallelogram`](@ref) with `AffineExpr` entries:
fully determined entries display as integers, undetermined entries
contain symbolic variables.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> P = hochschild_cohomology(Z);

julia> P[0, 0]  # h⁰(𝒪_Z) = 1
1

julia> P[3, 0]  # h⁰(∧³T_Z) = h⁰(ω_Z⁻¹) = 1 (CY3)
1
```
"""
function hochschild_cohomology(Z::ZeroLocus)
  X = Z.ambient
  E = Z.defining_bundle
  d = dimension(Z)
  var_counter = Ref(0)

  # ω_Z^{-1} lifted to X: by adjunction K_Z = (K_X ⊗ det E)|_Z,
  # so ω_Z^{-1} = (ω_X^{-1} ⊗ det(E)^{-1})|_Z.
  L_anti = tensor_product(anticanonical_bundle(X), dual(det(E)))
  L_can = tensor_product(canonical_bundle(X), det(E))

  is_fano = _is_fano_zero_locus(Z)

  E_dual = dual(E)

  # Precompute Sym^k(E*) and Ω^j_X
  syms = CompletelyReducibleBundle[symmetric_power(E_dual, k) for k in 0:d]
  omegas = CompletelyReducibleBundle[exterior_power(cotangent_bundle(X), j) for j in 0:d]

  # ── Build symbolic twisted Hodge matrices via LES ──────────────────────
  # M1[p+1, q+1] = h^q(Ω^p_Z ⊗ ω_Z^{-1}), symbolic AffineExpr
  # M2[p+1, q+1] = h^q(Ω^p_Z ⊗ ω_Z), symbolic AffineExpr
  M1 = Matrix{AffineExpr}(undef, d + 1, d + 1)
  M2 = Matrix{AffineExpr}(undef, d + 1, d + 1)
  for i in eachindex(M1)
    M1[i] = AffineExpr(0)
    M2[i] = AffineExpr(0)
  end

  for (M, L) in ((M1, L_anti), (M2, L_can))
    for p in 0:d
      if p == d
        # Ω^d_Z ⊗ L = (ω_X ⊗ det(E) ⊗ L)|_Z — direct Koszul
        F = tensor_product(tensor_product(canonical_bundle(X), det(E)), L)
        Hp = _restrict_to_zero_locus_les(Z, F, var_counter)
      else
        # Conormal filtration via symbolic LES chain
        # conormal_cohos[j+1] = H*(Z, Sym^{p-j}(E*) ⊗ Ω^j_X ⊗ L|_Z) for j=0..p
        conormal_cohos = Vector{AffineExpr}[]
        for j in 0:p
          F = tensor_product(tensor_product(syms[p - j + 1], omegas[j + 1]), L)
          Hj = _restrict_to_zero_locus_les(Z, F, var_counter)
          push!(conormal_cohos, Hj)
        end
        # Conormal chain: C_p = conormal_cohos[1] (highest Sym power)
        # SES: 0 → C_{j+1} → conormal_cohos[p-j+1] → C_j → 0
        # Reversed: terms[1] = K_p = conormal_cohos[1], terms[end] = K_0 = conormal_cohos[p+1]
        Hp = conormal_cohos[1]
        for k in 2:(p + 1)
          Hp = les_cokernel(Hp, conormal_cohos[k], var_counter)
        end
        # Apply vanishing: H^k(Z, F|_Z) = 0 for k > d
        if length(Hp) > d + 1
          mat = reshape(copy(Hp), 1, length(Hp))
          for k in (d + 1):(length(Hp) - 1)
            expr = mat[1, k + 1]
            is_zero_expr(expr) && continue
            _apply_equation_in_vars!(mat, expr)
          end
          Hp = AffineExpr[mat[1, k + 1] for k in 0:d]
        else
          Hp = Hp[1:(d + 1)]
        end
      end
      for q in 0:d
        M[p + 1, q + 1] = Hp[q + 1]
      end
    end
  end

  # ── Build symbolic polyvector data matrix ──────────────────────────────
  # data[p+1, q+1] = h^q(∧^p T_Z)
  #   = h^q(Ω^{d-p} ⊗ ω^{-1}) = M1[d-p+1, q+1]    (HKR identity)
  #   = h^{d-q}(Ω^p ⊗ ω)       = M2[p+1, d-q+1]    (Serre duality)
  data = Matrix{AffineExpr}(undef, d + 1, d + 1)
  for p in 0:d, q in 0:d
    data[p + 1, q + 1] = M1[d - p + 1, q + 1]
  end

  # ── Apply Serre duality cross-constraints ──────────────────────────────
  # data[p+1, q+1] = M2[p+1, d-q+1]
  constraint_changed = true
  while constraint_changed
    constraint_changed = false
    for p in 0:d, q in 0:d
      expr_hkr = data[p + 1, q + 1]   # = M1[d-p+1, q+1]
      expr_serre = M2[p + 1, d - q + 1]
      eq = expr_hkr - expr_serre
      if !isempty(eq.coeffs)
        constraint_changed = _apply_equation!(data, eq) || constraint_changed
        # Also propagate into M2
        _apply_equation!(M2, eq)
      end
    end

    # χ(∧^p T_Z) constraint: Σ_q (-1)^q data[p+1, q+1] = exact value
    for p in 0:d
      chi_p = _chi_omega_tensor(Z, d - p, L_anti)
      alt_sum = sum((-1)^q * data[p + 1, q + 1] for q in 0:d; init=AffineExpr(0))
      eq = alt_sum - AffineExpr(chi_p)
      if !isempty(eq.coeffs)
        constraint_changed = _apply_equation!(data, eq) || constraint_changed
      end
    end

    # Akizuki–Nakano vanishing for Fano: h^q(∧^p T) = 0 for q > p
    # Equivalently on twisted Hodge matrices:
    #   M1[p+1,q+1] = h^q(Ω^p ⊗ ω⁻¹) = 0 for p+q > d     (AN vanishing)
    #   M2[p+1,q+1] = h^q(Ω^p ⊗ ω)   = 0 for p+q < d     (Nakano dual)
    if is_fano
      for p in 0:d, q in (p + 1):d
        eq = data[p + 1, q + 1]
        if !is_zero_expr(eq)
          if !isempty(eq.coeffs)
            constraint_changed = _apply_equation!(data, eq) || constraint_changed
            _apply_equation!(M1, eq)
            _apply_equation!(M2, eq)
          end
          data[p + 1, q + 1] = AffineExpr(0)
          constraint_changed = true
        end
      end
      for p in 0:d, q in 0:d
        if p + q > d
          eq = M1[p + 1, q + 1]
          if !is_zero_expr(eq)
            if !isempty(eq.coeffs)
              constraint_changed = _apply_equation!(M1, eq) || constraint_changed
              _apply_equation!(data, eq)
              _apply_equation!(M2, eq)
            end
            M1[p + 1, q + 1] = AffineExpr(0)
            constraint_changed = true
          end
        end
        if p + q < d
          eq = M2[p + 1, q + 1]
          if !is_zero_expr(eq)
            if !isempty(eq.coeffs)
              constraint_changed = _apply_equation!(M2, eq) || constraint_changed
              _apply_equation!(data, eq)
              _apply_equation!(M1, eq)
            end
            M2[p + 1, q + 1] = AffineExpr(0)
            constraint_changed = true
          end
        end
      end
    end
  end

  # ── Fix boundary rows from exact Koszul computations ───────────────────
  # Row p=0: h^q(∧^0 T) = h^q(O_Z)
  (H_OZ, _) = cohomology_on_restriction(Z)
  for q in 0:d
    val = q <= H_OZ.dim_variety ? H_OZ[q] : BigInt(0)
    expr = data[1, q + 1]
    if !is_determined(expr) || expr.constant != val
      eq = expr - AffineExpr(val)
      !isempty(eq.coeffs) && _apply_equation!(data, eq)
      data[1, q + 1] = AffineExpr(val)
    end
  end
  # Row p=d: h^q(∧^d T) = h^q(ω_Z^{-1})
  (H_L, _) = cohomology_on_restriction(Z, L_anti)
  for q in 0:d
    val = q <= H_L.dim_variety ? H_L[q] : BigInt(0)
    expr = data[d + 1, q + 1]
    if !is_determined(expr) || expr.constant != val
      eq = expr - AffineExpr(val)
      !isempty(eq.coeffs) && _apply_equation!(data, eq)
      data[d + 1, q + 1] = AffineExpr(val)
    end
  end

  # ── Renumber remaining symbolic variables ──────────────────────────────
  _renumber_variables!(data)

  PolyvectorParallelogram(data, d)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge diamond printing
# ═══════════════════════════════════════════════════════════════════════════════

# Table format with no borders and no dividing lines, used for Hodge diamonds.
const _DIAMOND_FMT = TextTableFormat(;
  borders=text_table_borders__borderless,
  horizontal_line_at_beginning=false,
  horizontal_lines_at_column_labels=:none,
  horizontal_line_at_merged_column_labels=false,
  horizontal_line_after_column_labels=false,
  horizontal_lines_at_data_rows=:none,
  horizontal_line_before_row_group_label=false,
  horizontal_line_after_row_group_label=false,
  horizontal_line_after_data_rows=false,
  horizontal_line_before_summary_rows=false,
  horizontal_line_after_summary_rows=false,
  vertical_line_at_beginning=false,
  vertical_line_after_row_number_column=false,
  vertical_line_after_row_label_column=false,
  vertical_lines_at_data_columns=:none,
  vertical_line_after_data_columns=false,
  vertical_line_after_continuation_column=false,
)

"""
    print_hodge_diamond([io::IO,] h::Matrix)

Print a centred ASCII Hodge diamond for a variety of arbitrary dimension.

`h[p+1, q+1] = h^{p,q}` (1-based matrix indexing). The diamond is centred
using PrettyTables so that each column is as wide as its largest entry,
giving exact alignment regardless of digit count.

# Layout
`h^{p,q}` is placed at grid row `p+q+1`, column `p-q+d+1` in a
`(2d+1) × (2d+1)` string matrix (where `d = size(h,1) - 1`).

# Examples
```julia
julia> using PartialFlagVarieties

julia> X = Gr(2, 6);

julia> Z = zero_locus(reduce(direct_sum, [line_bundle(X, 1) for _ in 1:6]));

julia> print_hodge_diamond(stdout, hodge_numbers(Z));  # K3 of degree 14 in P^5
       1
    0      0
 1     20     1
    0      0
       1
```
"""
function print_hodge_diamond(io::IO, h::Matrix)
  d = size(h, 1) - 1
  sz = 2 * d + 1

  cells = fill("", sz, sz)
  for p in 0:d, q in 0:d
    cells[p + q + 1, p - q + d + 1] = string(h[p + 1, q + 1])
  end

  pretty_table(io, cells;
    table_format=_DIAMOND_FMT,
    show_column_labels=false,
    alignment=:c,
  )
end

print_hodge_diamond(h::Matrix) = print_hodge_diamond(stdout, h)
