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
  sum(((-1)^(p + q) * P[p, q] for p in 0:(P.dim), q in 0:(P.dim)); init=_pp_zero(T))
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

  TX = tangent_bundle(X)

  for p in 0:d
    wedge_p = exterior_power(TX, p)
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

For each row ``p \\le \\lfloor d/2 \\rfloor`` the conormal resolution of
``\\Omega^p_Z`` is restricted to ``Z`` through the Koszul resolution and the
long exact sequences are solved symbolically (see [`GradedConormal`](@ref)
and [`FilteredConormal`](@ref)).
On ambient varieties whose tangent bundle is not completely reducible, the
cotangent powers are treated as filtered bundles and their cohomology goes
through the spectral sequence of the filtration, so no silent degeneration
assumption is made.  Hodge symmetry, Serre duality, and Euler characteristic
constraints are then propagated to a fixed point, together with the
Lefschetz hyperplane theorem when the defining bundle splits off an ample
line bundle (see [`_lefschetz_inject!`](@ref)); the remaining rows follow
by Serre duality.

Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = \\mathrm{h}^{p,q}``.
Entries that are fully determined are plain integers (wrapped in `AffineExpr`);
entries that cannot be resolved from the available constraints contain
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
  d = dimension(Z)
  half = d ÷ 2
  var_counter = Ref(0)
  C = _conormal_data(Z, structure_sheaf(Z.ambient), half)

  M = fill(AffineExpr(0), d + 1, d + 1)
  for p in 0:half
    row = _conormal_row(C, p, var_counter)
    for q in 0:d
      M[p + 1, q + 1] = row[q + 1]
    end
  end

  chi_vals = BigInt[_chi_row(C, p) for p in 0:half]
  _propagate_hodge_constraints!(M, chi_vals, d)

  # The Lefschetz hyperplane theorem can pin down entries the long exact
  # sequences leave open.
  if _lefschetz_inject!(M, Z)
    _propagate_hodge_constraints!(M, chi_vals, d)
  end

  # Serre duality h^{p,q} = h^{d-p,d-q} fills the remaining rows.
  for p in (half + 1):d, q in 0:d
    M[p + 1, q + 1] = M[d - p + 1, d - q + 1]
  end

  _renumber_variables!(M)
end

"""
    hodge_numbers_symbolic(Z::ZeroLocus) -> Matrix{AffineExpr}

Symbolic version of `hodge_numbers`. When the long exact sequence does
not uniquely determine a Hodge number, the entry is an `AffineExpr`
involving symbolic variables ``x_0, x_1, \\ldots``.

This is a synonym of [`hodge_numbers`](@ref), which always returns the
symbolic form.

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
hodge_numbers_symbolic(Z::ZeroLocus) = hodge_numbers(Z)

"""
    hodge_numbers_les(Z::ZeroLocus) -> Matrix{AffineExpr}

Synonym of [`hodge_numbers`](@ref), kept for backwards compatibility.
"""
hodge_numbers_les(Z::ZeroLocus) = hodge_numbers(Z)

# ═══════════════════════════════════════════════════════════════════════════════
#  Twisted Hodge numbers of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
#  Conormal restriction of cotangent powers
#
#  Ω^p_Z is not the restriction of a bundle on X, so the Koszul resolution
#  does not apply to it directly.  Instead, the conormal sequence
#  0 → E^∨|_Z → Ω^1_X|_Z → Ω^1_Z → 0 yields the exact complex (exact because
#  the section is regular and Z smooth)
#
#      0 → Sym^p(E^∨) → Sym^{p-1}(E^∨) ⊗ Ω^1_X → ⋯ → Ω^p_X → Ω^p_Z → 0
#
#  on Z, whose remaining terms all are restrictions from X.  _conormal_row
#  computes H^*(Z, Ω^p_Z ⊗ L|_Z) by chaining this complex through
#  les_cokernel.
#
#  Restricting the terms Sym^{p-j}(E^∨) ⊗ Ω^j_X ⊗ L has two backends:
#
#   - GradedConormal: when the tangent bundle of X is completely reducible
#     (abelian nilradical), Ω^j_X equals its associated graded and everything
#     is done with deduplicated multiplicity dicts.
#   - FilteredConormal: otherwise Ω^j_X is a genuinely filtered bundle and
#     its cohomology goes through the spectral sequence of the filtration —
#     plain Borel–Weil–Bott on the associated graded would silently assume
#     degeneration.
#
#  The exact Euler characteristics χ(Ω^p_Z ⊗ L) only depend on K-theory
#  classes, so they always use the graded data.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    GradedConormal

Counts-level data for restricting the conormal terms
``\\mathrm{Sym}^{p-j}(E^\\vee) \\otimes \\Omega^j_X \\otimes L`` of a zero locus
``Z = Z(s) \\subset X`` to ``Z``, for ``p \\le p_{\\max}``, together with the
memoization state of the exact ``χ(\\Omega^p_Z \\otimes L)`` recursion.

Only valid as a restriction backend when the tangent bundle of ``X`` is
completely reducible; see [`FilteredConormal`](@ref) for the general case
(which still uses this type for χ and for the ``p = \\dim Z`` row).
"""
struct GradedConormal
  Z::ZeroLocus
  l_counts::Dict{IrrepLevi,Int}
  syms::Vector{Dict{IrrepLevi,Int}}    # gr Sym^k(E^∨), k = 0..pmax
  omegas::Vector{Dict{IrrepLevi,Int}}  # gr Ω^j_X,      j = 0..pmax
  wedge_counts::Vector{Dict{IrrepLevi,Int}}
  chi_memo::Dict{Tuple{Int,UInt},BigInt}
  chi_tp_memo::Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}
end

function GradedConormal(Z::ZeroLocus, L::CompletelyReducibleBundle, pmax::Int)
  X = Z.ambient
  E_dual = dual(Z.defining_bundle)
  GradedConormal(
    Z,
    _to_counts(L),
    Dict{IrrepLevi,Int}[_to_counts(symmetric_power(E_dual, k)) for k in 0:pmax],
    Dict{IrrepLevi,Int}[_to_counts(_cotangent_power(X, j)) for j in 0:pmax],
    Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)],
    Dict{Tuple{Int,UInt},BigInt}(),
    Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}(),
  )
end

"""
    FilteredConormal

Bundle-level backend for the conormal terms when ``\\Omega^1_X`` is genuinely
filtered (non-abelian nilradical): restriction to ``Z`` goes through the
spectral sequence of the filtration.  Wraps a [`GradedConormal`](@ref) for
everything that only depends on K-theory classes.
"""
struct FilteredConormal
  graded::GradedConormal
  syms_L::Vector{CompletelyReducibleBundle}  # Sym^k(E^∨) ⊗ L, k = 0..pmax
  omega_filt::Vector{FilteredBundle}         # filtered Ω^j_X, j = 0..jmax
end

_graded(C::GradedConormal) = C
_graded(C::FilteredConormal) = C.graded

"""
Choose the conormal restriction backend for rows ``p = 0, \\ldots, p_{\\max}``:
[`FilteredConormal`](@ref) when the tangent bundle of the ambient variety is
not completely reducible, [`GradedConormal`](@ref) otherwise.
"""
function _conormal_data(Z::ZeroLocus, L::CompletelyReducibleBundle, pmax::Int)
  graded = GradedConormal(Z, L, pmax)
  Omega_filt = dual(filtered_tangent_bundle(Z.ambient))
  n_filtration_steps(Omega_filt) > 1 || return graded

  # The row p = dim Z uses the ω_Z shortcut, so the filtered powers are only
  # needed for j ≤ min(pmax, dim Z - 1).
  jmax = min(pmax, dimension(Z) - 1)
  E_dual = dual(Z.defining_bundle)
  FilteredConormal(
    graded,
    CompletelyReducibleBundle[
      tensor_product(symmetric_power(E_dual, k), L) for k in 0:pmax
    ],
    FilteredBundle[exterior_power(Omega_filt, j) for j in 0:jmax],
  )
end

"""Restriction to ``Z`` of the conormal term ``\\mathrm{Sym}^{p-j}(E^\\vee) \\otimes \\Omega^j_X \\otimes L``."""
function _restrict_conormal_term(C::GradedConormal, p::Int, j::Int, var_counter::Ref{Int})
  f_counts = _tensor_product_counts(
    _tensor_product_counts(C.syms[p - j + 1], C.omegas[j + 1]), C.l_counts
  )
  _restrict_to_zero_locus_les(C.Z, f_counts, var_counter, C.wedge_counts)
end

function _restrict_conormal_term(C::FilteredConormal, p::Int, j::Int, var_counter::Ref{Int})
  term = tensor_product(C.omega_filt[j + 1], C.syms_L[p - j + 1])
  _restrict_to_zero_locus_les(_graded(C).Z, term, var_counter)
end

"""
    _conormal_row(C, p, var_counter) -> Vector{AffineExpr}

``H^q(Z, \\Omega^p_Z \\otimes L|_Z)`` for ``q = 0, \\ldots, \\dim Z``, where `C`
is a conormal backend from [`_conormal_data`](@ref).

Chains the conormal complex through `les_cokernel`, solving
``0 \\to C_{j-1} \\to \\mathrm{Sym}^{p-j}(E^\\vee) \\otimes \\Omega^j_X \\otimes L|_Z \\to C_j \\to 0``
for ``j = 1, \\ldots, p`` starting from ``C_0 = \\mathrm{Sym}^p(E^\\vee) \\otimes L|_Z``;
the last cokernel is ``\\Omega^p_Z \\otimes L|_Z``.  For ``p = \\dim Z`` the row is
the line bundle ``\\omega_Z \\otimes L`` and plain Koszul restriction is used.
"""
function _conormal_row(
  C::Union{GradedConormal,FilteredConormal}, p::Int, var_counter::Ref{Int}
)
  G = _graded(C)
  Z = G.Z
  d = dimension(Z)
  if p == d
    X = Z.ambient
    omega_Z_counts = _to_counts(tensor_product(canonical_bundle(X), det(Z.defining_bundle)))
    f_counts = _tensor_product_counts(omega_Z_counts, G.l_counts)
    return _restrict_to_zero_locus_les(Z, f_counts, var_counter, G.wedge_counts)
  end

  row = _restrict_conormal_term(C, p, 0, var_counter)
  for j in 1:p
    row = les_cokernel(row, _restrict_conormal_term(C, p, j, var_counter), var_counter)
  end

  # Vanishing H^k(Z, ·) = 0 for k > dim Z.
  _truncate_cohomology!(row, d)
end

"""Exact ``\\chi(Z, \\Omega^p_Z \\otimes L|_Z)`` from K-theory, memoized in `C`."""
function _chi_row(C::Union{GradedConormal,FilteredConormal}, p::Int)
  G = _graded(C)
  _chi_omega_tensor_counts_cached(
    G.Z, p, G.l_counts, G.chi_memo, G.omegas, G.wedge_counts, G.chi_tp_memo
  )
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Symbolic Hodge constraint helpers
# ═══════════════════════════════════════════════════════════════════════════════

_hodge_inconsistency() = error(
  "inconsistent Hodge constraints: is the section regular, with smooth zero locus " *
  "of the expected dimension?",
)

"""
Impose the exact value `M[i, j] = val`, eliminating a symbolic variable when
the entry is not determined yet.  A determined entry that disagrees signals
inconsistent input data and raises an error.
"""
function _inject_exact!(M::Matrix{AffineExpr}, i::Int, j::Int, val::BigInt)
  expr = M[i, j]
  if is_determined(expr)
    expr.constant == val || _hodge_inconsistency()
    return false
  end
  _apply_equation!(M, expr - AffineExpr(val))
  M[i, j] = AffineExpr(val)
  true
end

"""Impose `M[pi, qi] = target`; see [`_inject_exact!`](@ref)."""
function _apply_hodge_constraint!(M::Matrix{AffineExpr}, pi::Int, qi::Int, target::BigInt)
  expr = M[pi, qi]
  if is_determined(expr)
    expr.constant == target || _hodge_inconsistency()
    return false
  end
  _apply_equation!(M, expr - AffineExpr(target))
end

"""Impose `M[p1, q1] = M[p2, q2]`, erroring on a determined contradiction."""
function _apply_hodge_pair!(M::Matrix{AffineExpr}, p1::Int, q1::Int, p2::Int, q2::Int)
  eq = M[p1, q1] - M[p2, q2]
  if is_determined(eq)
    eq.constant == 0 || _hodge_inconsistency()
    return false
  end
  _apply_equation!(M, eq)
end

"""
Propagate Hodge symmetry, Serre duality, and Euler characteristic constraints
through the ``(d+1) \\times (d+1)`` matrix `M` of ``h^{p,q}`` candidates, of
which only rows ``p = 0, \\ldots, \\lfloor d/2 \\rfloor`` are filled (the others
follow by Serre duality afterwards).  `chi_vals[p+1]` is the exact
``\\chi(\\Omega^p_Z)``.  Iterates to a fixed point.
"""
function _propagate_hodge_constraints!(
  M::Matrix{AffineExpr}, chi_vals::Vector{BigInt}, d::Int
)
  half = d ÷ 2
  changed = true
  while changed
    changed = false

    # Hodge symmetry against the exact row 0: h^{p,0} = h^{0,p}, and
    # Serre duality h^{p,d} = h^{d-p,0} = h^{0,d-p}.
    for p in 1:half
      if is_determined(M[1, p + 1])
        changed |= _apply_hodge_constraint!(M, p + 1, 1, M[1, p + 1].constant)
      end
      dp = d - p
      if dp != p && is_determined(M[1, dp + 1])
        changed |= _apply_hodge_constraint!(M, p + 1, d + 1, M[1, dp + 1].constant)
      end
    end

    # χ(Ω^p_Z) per computed row.
    for p in 0:half
      eq = _alternating_sum(M, p + 1, d) - AffineExpr(chi_vals[p + 1])
      is_determined(eq) && eq.constant != 0 && _hodge_inconsistency()
      changed |= _apply_equation!(M, eq)
    end

    # Hodge symmetry h^{p,q} = h^{q,p} within the computed rows.
    for p in 0:half, q in 0:half
      p == q && continue
      changed |= _apply_hodge_pair!(M, p + 1, q + 1, q + 1, p + 1)
    end

    # Middle-row Serre duality h^{half,q} = h^{half,d-q} for even d.
    if d % 2 == 0
      for q in 0:(half - 1)
        changed |= _apply_hodge_pair!(M, half + 1, q + 1, half + 1, d - q + 1)
      end
    end

    # Hodge combined with Serre: h^{p,q} = h^{d-q,d-p} when both rows are
    # among the computed ones.
    for p in 0:half, q in 0:d
      dq = d - q
      0 <= dq <= half || continue
      (p == dq && q == d - p) && continue
      changed |= _apply_hodge_pair!(M, p + 1, q + 1, dq + 1, d - p + 1)
    end
  end
  M
end

"""
    _lefschetz_inject!(M, Z) -> Bool

Resolve entries of the Hodge matrix via the Lefschetz hyperplane theorem.

When ``E = E' \\oplus L`` with ``L`` an ample line bundle, ``Z = Z(E)`` is an
ample divisor in ``Z' = Z(E')``, so ``h^{p,q}(Z) = h^{p,q}(Z')`` for
``p + q < \\dim Z``.  The Hodge numbers of ``Z'`` are computed recursively
(further ample line-bundle summands strip off the same way; with no summand
left they are the Hodge numbers of the ambient space) and the determined
values are injected into the computed rows ``p \\le \\lfloor d/2 \\rfloor``.

Only attempted when undetermined entries remain below the middle total
degree, since the recursive computation is not free.  Returns whether
anything was injected.
"""
function _lefschetz_inject!(M::Matrix{AffineExpr}, Z::ZeroLocus)
  d = dimension(Z)
  half = d ÷ 2
  any(
    !is_determined(M[p + 1, q + 1]) for p in 0:half for q in 0:d if p + q < d
  ) || return false

  X = Z.ambient
  comps = components(Z.defining_bundle)
  i = findfirst(
    c -> fiber_dimension(c) == 1 &&
      is_ample_line_bundle(CompletelyReducibleBundle(X, [c])),
    comps,
  )
  i === nothing && return false

  rest = comps[setdiff(eachindex(comps), i)]
  H = if isempty(rest)
    hodge_numbers(X)
  else
    hodge_numbers(zero_locus(CompletelyReducibleBundle(X, rest)))
  end

  changed = false
  for p in 0:half, q in 0:d
    p + q < d || continue
    v = H[p + 1, q + 1]
    v isa AffineExpr && !is_determined(v) && continue
    val = v isa AffineExpr ? v.constant : BigInt(v)
    changed |= _inject_exact!(M, p + 1, q + 1, val)
  end
  changed
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Twisted Hodge engine
# ═══════════════════════════════════════════════════════════════════════════════

function _twisted_hodge_symbolic(
  Z::ZeroLocus, L::CompletelyReducibleBundle, var_counter::Ref{Int}
)
  X = Z.ambient
  E = Z.defining_bundle
  d = dimension(Z)
  C = _conormal_data(Z, L, d)

  M = fill(AffineExpr(0), d + 1, d + 1)
  for p in 0:d
    row = _conormal_row(C, p, var_counter)
    for q in 0:d
      M[p + 1, q + 1] = row[q + 1]
    end
  end

  # Exact boundary rows, when plain Koszul restriction determines them:
  # p = 0 is L|_Z and p = d is (ω_Z ⊗ L)|_Z.
  (H_L, det_L) = cohomology_on_restriction(Z, L)
  if det_L
    for q in 0:d
      _inject_exact!(M, 1, q + 1, BigInt(H_L[q]))
    end
  end
  omega_Z_L = tensor_product(tensor_product(canonical_bundle(X), det(E)), L)
  (H_top, det_top) = cohomology_on_restriction(Z, omega_Z_L)
  if det_top
    for q in 0:d
      _inject_exact!(M, d + 1, q + 1, BigInt(H_top[q]))
    end
  end

  # Kodaira–Akizuki–Nakano vanishing: for ample L we get
  # H^q(Z, Ω^p_Z ⊗ L|_Z) = 0 for p + q > d, and for antiample L the Serre
  # dual statement gives vanishing for p + q < d.  Ampleness of L on X
  # restricts to ampleness on Z.
  if is_ample_line_bundle(L)
    for p in 0:d, q in (d - p + 1):d
      _inject_exact!(M, p + 1, q + 1, BigInt(0))
    end
  elseif is_ample_line_bundle(dual(L))
    for p in 0:d, q in 0:(d - p - 1)
      _inject_exact!(M, p + 1, q + 1, BigInt(0))
    end
  end

  # χ constraint per row: Σ_q (-1)^q h^q(Ω^p_Z ⊗ L) is exact from K-theory.
  for p in 0:d
    eq = _alternating_sum(M, p + 1, d) - AffineExpr(_chi_row(C, p))
    is_determined(eq) && eq.constant != 0 && _hodge_inconsistency()
    _apply_equation!(M, eq)
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
  M = _twisted_hodge_symbolic(Z, L, Ref(0))
  _renumber_variables!(M)
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
  # Two families of constraints, iterated to a fixed point (Akizuki–Nakano
  # vanishing for confirmed Fano Z is already applied inside the twisted
  # engine, since L_anti is then an ample line bundle):
  #   1. Serre duality:  data[p+1, q+1] = M2[p+1, d-q+1]
  #   2. Euler char:     Σ_q (-1)^q h^q(∧^p T) = χ(Ω^{d-p} ⊗ ω⁻¹)
  # χ(∧^p T_Z) = χ(Ω^{d-p}_Z ⊗ ω_Z^{-1}), exact from the conormal recursion.
  l_anti_counts = _to_counts(L_anti)
  chi_memo = Dict{Tuple{Int,UInt},BigInt}()
  chi_tp_memo = Dict{Tuple{UInt,Int,Bool},Dict{IrrepLevi,Int}}()
  chi_omegas = Dict{IrrepLevi,Int}[_to_counts(_cotangent_power(X, j)) for j in 0:d]
  wedge_counts = Dict{IrrepLevi,Int}[_to_counts(w) for w in _koszul_wedges!(Z)]
  chi_polyvector = BigInt[
    _chi_omega_tensor_counts_cached(
      Z, d - p, l_anti_counts, chi_memo, chi_omegas, wedge_counts, chi_tp_memo
    ) for p in 0:d
  ]

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
    for p in 0:d
      eq = _alternating_sum(data, p + 1, d) - AffineExpr(chi_polyvector[p + 1])
      if !isempty(eq.coeffs)
        constraint_changed = _apply_equation!(data, eq) || constraint_changed
      end
    end

    # Akizuki–Nakano vanishing (for confirmed Fano) is already applied
    # inside _twisted_hodge_symbolic: L_anti is then an ample line bundle, so
    # M1 vanishes for p + q > d (equivalently h^q(∧^p T_Z) = 0 for q > p)
    # and M2, twisted by the antiample ω_Z, vanishes for p + q < d.
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
