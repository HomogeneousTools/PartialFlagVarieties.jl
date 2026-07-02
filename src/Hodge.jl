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
export hochschild_cohomology
export PolyvectorParallelogram
export print_hodge_diamond

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge numbers  h^{p,q}(X)  for partial flag varieties
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers(X::PartialFlagVariety) -> Matrix{BigInt}

Compute the Hodge numbers ``\\mathrm{h}^{p,q}(X) = \\dim H^q(X, \\Omega^p_X)``
of the partial flag variety ``X = G/P``.

Since ``G/P`` is rational (and simply connected), the Hodge diamond is
diagonal: ``\\mathrm{h}^{p,q} = 0`` for ``p \\neq q``, and ``\\mathrm{h}^{p,p} = b_{2p}``
where ``b_i`` are the Betti numbers.

Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = \\mathrm{h}^{p,q}``,
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

Here `j` means tensoring by the generator of ``\\operatorname{Pic}(X)``. This
therefore requires ``X`` to have Picard rank 1, equivalently to be a
generalized Grassmannian.

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
    PolyvectorParallelogram{T}

The Hochschild–Kostant–Rosenberg decomposition of Hochschild cohomology:

``\\mathrm{HH}^n(X) = \\bigoplus_{p+q=n} H^q(X, \\bigwedge^p T_X)``

Stored as a matrix where entry ``[p+1, q+1] = h^q(X, \\bigwedge^p T_X)``.
The type parameter `T` is `BigInt` when all entries are determined and
[`AffineExpr`](@ref) when some entries depend on undetermined connecting-map
ranks from long exact sequences (typically for zero loci).

Access the mathematically labelled entry with `P[p, q]`, which uses 0-based
polyvector degree `p` and cohomological degree `q`.

# Fields
- `data::Matrix{T}`: the HKR decomposition matrix
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
function Base.getindex(P::PolyvectorParallelogram{T}, p::Int, q::Int) where {T}
  0 <= p <= P.dim || return _pp_zero(T)
  0 <= q <= P.dim || return _pp_zero(T)
  P.data[p + 1, q + 1]
end

"""
    euler_characteristic(P::PolyvectorParallelogram)

Compute the Euler characteristic of Hochschild cohomology:
``\\chi(\\mathrm{HH}^*) = \\sum_{p,q} (-1)^{p+q} h^q(X, \\bigwedge^p T_X)``
"""
function euler_characteristic(P::PolyvectorParallelogram{T}) where {T}
  result = _pp_zero(T)
  for p in 0:(P.dim)
    for q in 0:(P.dim)
      result += (-1)^(p + q) * P[p, q]
    end
  end
  result
end

"""
    getindex(P::PolyvectorParallelogram, n::Int)

Return ``\\dim \\mathrm{HH}^n(X) = \\sum_{p+q=n} h^q(X, \\bigwedge^p T_X)``.
"""
function Base.getindex(P::PolyvectorParallelogram{T}, n::Int) where {T}
  result = _pp_zero(T)
  for p in 0:(P.dim)
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
Use `P[p, q]` to read the entry ``h^q(X, \\bigwedge^p T_X)``.

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

Compute the Hodge diamond ``\\mathrm{h}^{p,q}(Z)`` for ``p, q = 0, \\ldots, \\dim Z``.

Uses the Koszul resolution and long exact sequences.
Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = \\mathrm{h}^{p,q}``.
Entries that are fully determined are plain integers (wrapped in `AffineExpr`);
entries that cannot be resolved from Koszul + symmetry constraints contain
symbolic variables.

As throughout the zero-locus API, this assumes that `Z` is the smooth zero
locus of a regular section.

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

"""
    hodge_numbers_symbolic(Z::ZeroLocus) -> Matrix{AffineExpr}

Symbolic version of `hodge_numbers`. When the long exact sequence does
not uniquely determine a Hodge number, the entry is an `AffineExpr`
involving symbolic variables ``x_0, x_1, \\ldots``.

This function delegates to [`hodge_numbers_les`](@ref), which is the
shared symbolic Hodge solver used for zero loci.

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
  hodge_numbers_les(Z)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Twisted Hodge numbers of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

# ── Shared symbolic engine for twisted Hodge computation ─────────────────────
#
# _twisted_hodge_symbolic(Z, L, var_counter) computes h^q(Ω^p_Z ⊗ L|_Z)
# symbolically via:
#   1. Koszul resolution for each conormal graded piece (tries numeric first)
#   2. les_cokernel to chain the conormal SES
#   3. χ(Ω^p_Z ⊗ L) constraint per row (exact from _chi_omega_tensor_counts)
#   4. Exact boundary injection for p = 0 and p = d
#
# Both twisted_hodge_numbers and hochschild_cohomology delegate to this.

function _twisted_hodge_symbolic(
  Z::ZeroLocus, L::CompletelyReducibleBundle, var_counter::Ref{Int}
)
  X = Z.ambient
  E = Z.defining_bundle
  E_dual = dual(E)
  d = dimension(Z)

  # Precompute Sym^k(E*) and Ω^j_X as multiplicity dicts for k,j = 0..d.
  # The _to_counts representation (Dict{IrrepLevi,Int}) deduplicates
  # identical summands, avoiding redundant tensor-product calls.
  syms_counts = Dict{IrrepLevi,Int}[
    _to_counts(symmetric_power(E_dual, k)) for k in 0:d
  ]
  omegas_counts = Dict{IrrepLevi,Int}[
    _to_counts(_cotangent_power(X, j)) for j in 0:d
  ]
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]

  # ω_Z lifted to X (before restriction), as a multiplicity dict for the
  # p = d branch (Ω^d_Z = ω_Z).
  omega_Z_counts = _to_counts(tensor_product(canonical_bundle(X), det(E)))

  # Convert the twist line bundle to a multiplicity dict for
  # _tensor_product_counts.
  l_counts = _to_counts(L)

  M = fill(AffineExpr(0), d + 1, d + 1)

  for p in 0:d
    if p == d
      # Ω^d_Z ⊗ L = (ω_Z ⊗ L)|_Z — use the precomputed ω_Z counts.
      f_counts = _tensor_product_counts(omega_Z_counts, l_counts)
      Hp = _restrict_to_zero_locus_les(Z, f_counts, var_counter, wedge_counts)
    else
      # ── Conormal filtration ─────────────────────────────────────────
      # For 0 < p < d, Ω^p_Z is intrinsic to Z — it is NOT the
      # restriction of a bundle from X, so we cannot apply the Koszul
      # resolution directly.  Instead, the conormal exact sequence
      #   0 → E∨|_Z → Ω¹_X|_Z → Ω¹_Z → 0
      # induces a filtration on Ω^p_Z whose graded pieces ARE
      # restrictions from X:
      #   Gr_j = Sym^{p-j}(E∨) ⊗ Ω^j_X   for j = 0, …, p.
      #
      # We compute H*(Z, Gr_j ⊗ L|_Z) for each graded piece via
      # the Koszul resolution, then recover H*(Z, Ω^p_Z ⊗ L|_Z)
      # by chaining short exact sequences:
      #   0 → C_{k-1} → H*(Z, Gr_k ⊗ L|_Z) → C_k → 0
      # with C_0 = H*(Z, Sym^p(E∨) ⊗ L|_Z).  After p steps,
      # C_p = H*(Z, Ω^p_Z ⊗ L|_Z) is the desired cohomology.
      conormal_cohomologies = Vector{AffineExpr}[]
      for j in 0:p
        f_counts = _tensor_product_counts(
          _tensor_product_counts(syms_counts[p - j + 1], omegas_counts[j + 1]),
          l_counts,
        )
        Hj = _restrict_to_zero_locus_les(Z, f_counts, var_counter, wedge_counts)
        push!(conormal_cohomologies, Hj)
      end
      # les_cokernel(A, B) solves 0 → A → B → C → 0 for H*(C).
      Hp = conormal_cohomologies[1]
      for k in 2:(p + 1)
        Hp = les_cokernel(Hp, conormal_cohomologies[k], var_counter)
      end
      # Enforce vanishing H^k(Z, F|_Z) = 0 for k > dim Z.
      Hp = copy(Hp)
      for k in (d + 1):(length(Hp) - 1)
        is_zero_expr(Hp[k + 1]) || _apply_equation!(Hp, Hp[k + 1])
      end
      Hp = Hp[1:(d + 1)]
    end
    for q in 0:d
      M[p + 1, q + 1] = Hp[q + 1]
    end
  end

  # ── Inject exact boundary values for p = 0 and p = d ─────────────────
  # p = 0: h^q(Ω⁰_Z ⊗ L) = h^q(L|_Z), computable exactly.
  (H_L, _) = cohomology_on_restriction(Z, L)
  for q in 0:d
    val = BigInt(q <= H_L.dim_variety ? H_L[q] : 0)
    _inject_exact!(M, 1, q + 1, val)
  end
  # p = d: h^q(Ω^d_Z ⊗ L) = h^q((ω_Z ⊗ L)|_Z), computable exactly.
  omega_Z_L = tensor_product(tensor_product(canonical_bundle(X), det(E)), L)
  (H_top, _) = cohomology_on_restriction(Z, omega_Z_L)
  for q in 0:d
    val = BigInt(q <= H_top.dim_variety ? H_top[q] : 0)
    _inject_exact!(M, d + 1, q + 1, val)
  end

  # ── χ constraint per row ─────────────────────────────────────────────
  # Σ_q (-1)^q h^q(Ω^p_Z ⊗ L) = χ(Ω^p_Z ⊗ L), exact from the Koszul complex.
  chi_memo = Dict{Tuple{Int,UInt},BigInt}()
  chi_tp_memo = Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}()
  for p in 0:d
    chi_p = _chi_omega_tensor_counts_cached(
      Z, p, l_counts, chi_memo, omegas_counts, wedge_counts, chi_tp_memo
    )
    alt_sum = _alternating_sum(M, p + 1, d)
    eq = alt_sum - AffineExpr(chi_p)
    !isempty(eq.coeffs) && _apply_equation!(M, eq)
  end

  M
end

"""
    twisted_hodge_numbers(Z::ZeroLocus, L::CompletelyReducibleBundle) -> Matrix{AffineExpr}

Compute the twisted Hodge numbers ``h^q(Z, \\Omega^p_Z \\otimes L|_Z)``
for the zero locus ``Z = Z(s) \\subset X = G/P`` and a line bundle ``L``
on ``X``.

Uses the conormal sequence to decompose ``\\Omega^p_Z`` into restrictions
from the ambient variety, computes each piece via the Koszul resolution,
and chains the results symbolically.

Returns a ``(d+1) \\times (d+1)`` matrix of `AffineExpr` entries where
``d = \\dim Z`` and entry ``[p+1, q+1] = h^q(Z, \\Omega^p_Z \\otimes L|_Z)``.
Fully determined entries display as integers; use `is_determined(M[p+1,q+1])`
to check.

The bundle `L` is given on the ambient variety; the computation is for its
restriction to `Z`.

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
  _twisted_hodge_symbolic(Z, L, Ref(0))
end

"""
    twisted_hodge_numbers(Z::ZeroLocus, j::Integer) -> Matrix{AffineExpr}

Compute the twisted Hodge numbers ``h^q(Z, \\Omega^p_Z(j))`` where the
twist is by ``\\mathcal{O}_X(j)|_Z``.

Here `j` denotes the integer power of the Picard-rank-1 generator on the
ambient variety. This therefore requires the ambient variety to have Picard
rank 1.

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

"""
Replace the symbolic AffineExpr in matrix M at (i, j) with the exact
integer value `val`, propagating the resulting equation into `M`.
"""
function _inject_exact!(M::Matrix{AffineExpr}, i::Int, j::Int, val::BigInt)
  expr = M[i, j]
  if !is_determined(expr) || expr.constant != val
    eq = expr - AffineExpr(val)
    !isempty(eq.coeffs) && _apply_equation!(M, eq)
    M[i, j] = AffineExpr(val)
  end
end

"""
    hochschild_cohomology(Z::ZeroLocus) -> PolyvectorParallelogram{AffineExpr}

Compute the Hochschild cohomology of ``Z`` via the HKR decomposition:

``\\mathrm{HH}^n(Z) = \\bigoplus_{p+q=n} H^q(Z, \\bigwedge^p T_Z)``

Uses the identity ``H^q(Z, \\bigwedge^p T_Z) = H^q(Z, \\Omega^{d-p}_Z
\\otimes \\omega_Z^{-1})`` where ``d = \\dim Z``, reducing the computation
to twisted Hodge numbers with the anticanonical twist
``L = \\omega_X^{-1} \\otimes \\det(E)^{-1}`` lifted to the ambient variety.

Akizuki–Nakano vanishing (``h^q(\\bigwedge^p T_Z) = 0`` for ``q > p``) is
applied when ``\\omega_Z^{-1}`` is *confirmed* ample from the ambient ``G/P``
(i.e. all Picard-coordinate differences are strictly positive).  When some
coordinate is zero the Fano status of ``Z`` is undetermined and vanishing is
not assumed.

Returns a [`PolyvectorParallelogram`](@ref) with `AffineExpr` entries:
fully determined entries display as integers, undetermined entries
contain symbolic variables.

As for the homogeneous case, use `P[p, q]` for the entry
``h^q(Z, \\bigwedge^p T_Z)``.

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

  # ── Lift ω_Z and ω_Z⁻¹ to the ambient variety ─────────────────────────
  # By adjunction, K_Z = (K_X ⊗ det E)|_Z.  We need two line bundles on X
  # whose restrictions to Z give ω_Z⁻¹ and ω_Z:
  #   L_anti = ω_X⁻¹ ⊗ det(E)⁻¹  →  restricts to ω_Z⁻¹
  #   L_can  = ω_X  ⊗ det(E)      →  restricts to ω_Z
  # These are dual: L_can = L_anti∨.
  L_anti = tensor_product(anticanonical_bundle(X), dual(det(E)))
  L_can = dual(L_anti)

  is_fano = is_strongly_fano(Z)

  # ── Build two symbolic twisted Hodge matrices ──────────────────────────
  # We compute h^q(Ω^p_Z ⊗ L) for L = ω_Z⁻¹ and L = ω_Z independently,
  # sharing a single var_counter so all symbolic variables are distinct.
  # Each call to _twisted_hodge_symbolic already applies: conormal
  # filtration, boundary injection (p=0 and p=d), and χ constraints.
  #
  #   M1[p+1, q+1] = h^q(Ω^p_Z ⊗ ω_Z⁻¹)   (anticanonical twist)
  #   M2[p+1, q+1] = h^q(Ω^p_Z ⊗ ω_Z)      (canonical twist)
  #
  # Via the HKR isomorphism ∧^p T_Z ≅ Ω^{d-p}_Z ⊗ ω_Z⁻¹:
  #   h^q(∧^p T_Z) = M1[d-p+1, q+1]
  #
  # Via Serre duality h^q(∧^p T_Z) = h^{d-q}((∧^p T_Z)∨ ⊗ ω_Z)
  #                                 = h^{d-q}(Ω^p_Z ⊗ ω_Z):
  #   h^q(∧^p T_Z) = M2[p+1, d-q+1]
  #
  # M1 and M2 have different symbolic variables (from the shared counter).
  # Cross-constraining these two representations resolves entries that
  # neither matrix alone can determine.
  M1 = _twisted_hodge_symbolic(Z, L_anti, var_counter)
  M2 = _twisted_hodge_symbolic(Z, L_can, var_counter)

  # ── Additional boundary injection for M2 via Serre duality ─────────────
  # M2[1, q+1] = h^q(Ω⁰_Z ⊗ ω_Z) = h^q(ω_Z) = h^{d-q}(𝒪_Z), and the
  # h^q(𝒪_Z) values are already exact in M1[d+1, :] (its p=d row with
  # L=ω_Z⁻¹ gives h^q(Ω^d ⊗ ω_Z⁻¹) = h^q(𝒪_Z)).
  for q in 0:d
    val_OZ = M1[d + 1, d - q + 1]
    if is_determined(val_OZ)
      _inject_exact!(M2, 1, q + 1, BigInt(val_OZ.constant))
    end
  end

  # ── Build symbolic polyvector data matrix ──────────────────────────────
  # h^q(∧^p T_Z) = h^q(Ω^{d-p}_Z ⊗ ω_Z⁻¹) = M1[d-p+1, q+1]  (HKR)
  data = Matrix{AffineExpr}(undef, d + 1, d + 1)
  for p in 0:d, q in 0:d
    data[p + 1, q + 1] = M1[d - p + 1, q + 1]
  end

  # ── Constraint propagation loop ────────────────────────────────────────
  # Three families of constraints, iterated to a fixed point:
  #   1. Serre duality:  data[p+1, q+1] = M2[p+1, d-q+1]
  #   2. Euler char:     Σ_q (-1)^q h^q(∧^p T) = χ(Ω^{d-p} ⊗ ω⁻¹)
  #   3. Akizuki–Nakano: h^q(∧^p T) = 0 for q > p  (Fano only)
  l_anti_counts = _to_counts(L_anti)
  constraint_changed = true
  while constraint_changed
    constraint_changed = false

    # 1. Serre duality: h^q(∧^p T) = h^{d-q}(Ω^p ⊗ ω)
    for p in 0:d, q in 0:d
      expr_hkr = data[p + 1, q + 1]
      expr_serre = M2[p + 1, d - q + 1]
      eq = expr_hkr - expr_serre
      if !isempty(eq.coeffs)
        constraint_changed = _apply_equation!(data, eq) || constraint_changed
        _apply_equation!(M2, eq)
      end
    end

    # 2. χ(∧^p T_Z) = exact value from conormal recursion
    chi_memo = Dict{Tuple{Int,UInt},BigInt}()
    chi_tp_memo = Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}()
    chi_omegas_counts = Dict{IrrepLevi,Int}[
      _to_counts(_cotangent_power(X, j)) for j in 0:d
    ]
    wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
    for p in 0:d
      chi_p = _chi_omega_tensor_counts_cached(
        Z, d - p, l_anti_counts, chi_memo, chi_omegas_counts, wedge_counts, chi_tp_memo
      )
      alt_sum = _alternating_sum(data, p + 1, d)
      eq = alt_sum - AffineExpr(chi_p)
      if !isempty(eq.coeffs)
        constraint_changed = _apply_equation!(data, eq) || constraint_changed
      end
    end

    # 3. Akizuki–Nakano vanishing for Fano zero loci
    #   h^q(∧^p T) = 0 for q > p  (on data)
    #   h^q(Ω^p ⊗ ω⁻¹) = 0 for p + q > d  (on M1)
    #   h^q(Ω^p ⊗ ω)   = 0 for p + q < d  (Nakano dual, on M2)
    if is_fano === true
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

`h[p+1, q+1] = \\mathrm{h}^{p,q}` (1-based matrix indexing). The diamond is centred
using PrettyTables so that each column is as wide as its largest entry,
giving exact alignment regardless of digit count.

# Layout
`\\mathrm{h}^{p,q}` is placed at grid row `p+q+1`, column `p-q+d+1` in a
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
