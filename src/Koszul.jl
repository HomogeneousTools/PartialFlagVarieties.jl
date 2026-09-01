# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul.jl — Long exact sequence algebra for Koszul resolutions
#
#  Pure algebra module: solves long exact sequences arising from
#  Koszul-type filtrations.  Knows nothing about flag varieties or bundles.
#
#  Given the dimension-valued cohomology of two terms of a short exact
#  sequence, computes the cohomology of the third term as far as exactness
#  determines it.  Undetermined connecting-map ranks are either reported
#  (determined-dimension solvers, which return a flag) or turned into
#  symbolic variables (the `AffineExpr`-valued solvers).
#
#  Layout:
#   1. AffineExpr        — integer affine expressions c + Σ kⱼ xⱼ, with
#                          substitution, elimination, and renumbering.
#   2. Intervals         — bound-consistency propagation over systems of
#                          nonnegative affine quantities.
#   3. BigInt solvers    — connecting-rank bound propagation with a
#                          `determined` flag.
#   4. δ-based solvers   — one symbolic variable per open connecting rank
#                          (solve_ses_cohomology_symbolic and the filtration
#                          chain).
#   5. Entry-based solvers — one symbolic variable per unknown entry, with
#                          equations *and inequalities* harvested from
#                          exactness (les_cokernel, les_kernel).
#
#  Conventions: short exact sequences are written 0 → A → B → C → 0, with
#  aᵢ = dim Hⁱ(A) and so on; δᵢ = rank(Hⁱ(C) → Hⁱ⁺¹(A)) is the connecting
#  rank, so exactness reads cᵢ = bᵢ - aᵢ + δᵢ₋₁ + δᵢ.  Every symbolic
#  variable is a nonnegative integer (a dimension, a rank, or a rank shifted
#  by a known lower bound).
# ═══════════════════════════════════════════════════════════════════════════════

export solve_ses_cohomology, solve_koszul_filtration
export AffineExpr, is_determined, is_zero_expr, symbolic_variable
export solve_ses_cohomology_symbolic, solve_koszul_filtration_symbolic

# ═══════════════════════════════════════════════════════════════════════════════
#  AffineExpr — affine expressions for symbolic cohomology
# ═══════════════════════════════════════════════════════════════════════════════

"""
    AffineExpr

An affine expression ``c + \\sum_j k_j x_j`` representing a cohomology
dimension that depends on undetermined connecting-map ranks.

Each ``x_j`` corresponds to a free parameter in a long exact sequence
whose rank could not be uniquely determined.
"""
struct AffineExpr
  constant::BigInt
  coeffs::Dict{Int,BigInt}  # var_id => coefficient; only non-zero entries
end

AffineExpr(c::Integer) = AffineExpr(BigInt(c), Dict{Int,BigInt}())

"""Create a fresh symbolic variable ``x_i``."""
symbolic_variable(var_id::Int) = AffineExpr(
  BigInt(0), Dict{Int,BigInt}(var_id => BigInt(1))
)

"""
Create a fresh symbolic variable and advance the counter.  Ids only need to
be distinct within one computation; the user-facing output is compacted by
`_renumber_variables!` at the end.
"""
function _fresh_variable(var_counter::Ref{Int})
  x = symbolic_variable(var_counter[])
  var_counter[] += 1
  x
end

"""Check whether the expression is fully determined (no symbolic variables)."""
is_determined(e::AffineExpr) = isempty(e.coeffs)

"""Return whether every entry of a cohomology object is determined."""
is_determined(::Cohomology) = true
is_determined(H::Cohomology{AffineExpr}) = all(is_determined, H.entries)

"""Check whether the expression is identically zero."""
is_zero_expr(e::AffineExpr) = e.constant == 0 && isempty(e.coeffs)

_determined_bigints(entries::Vector{AffineExpr}) = BigInt[e.constant for e in entries]

_as_affine(v::Vector{AffineExpr}) = copy(v)
_as_affine(v::Vector{BigInt}) = AffineExpr.(v)

# ─── Arithmetic ──────────────────────────────────────────────────────────────
#
# AffineExpr is treated as immutable throughout: every operation returns a
# fresh expression (sharing is safe, mutation never happens in place).

function _merge_coeffs(op, a::Dict{Int,BigInt}, b::Dict{Int,BigInt})
  merged = copy(a)
  for (var, coeff) in b
    merged[var] = op(get(merged, var, BigInt(0)), coeff)
    merged[var] == 0 && delete!(merged, var)
  end
  merged
end

function Base.:+(a::AffineExpr, b::AffineExpr)
  constant = a.constant + b.constant
  isempty(a.coeffs) && isempty(b.coeffs) && return AffineExpr(constant)
  isempty(a.coeffs) && return AffineExpr(constant, copy(b.coeffs))
  isempty(b.coeffs) && return AffineExpr(constant, copy(a.coeffs))
  AffineExpr(constant, _merge_coeffs(+, a.coeffs, b.coeffs))
end

function Base.:-(a::AffineExpr, b::AffineExpr)
  constant = a.constant - b.constant
  isempty(a.coeffs) && isempty(b.coeffs) && return AffineExpr(constant)
  isempty(b.coeffs) && return AffineExpr(constant, copy(a.coeffs))
  isempty(a.coeffs) &&
    return AffineExpr(constant, Dict{Int,BigInt}(k => -v for (k, v) in b.coeffs))
  AffineExpr(constant, _merge_coeffs(-, a.coeffs, b.coeffs))
end

function Base.:-(a::AffineExpr)
  isempty(a.coeffs) && return AffineExpr(-a.constant)
  AffineExpr(-a.constant, Dict{Int,BigInt}(k => -v for (k, v) in a.coeffs))
end

function Base.:*(c::Integer, a::AffineExpr)
  c == 0 && return AffineExpr(0)
  c == 1 && return a
  c == -1 && return -a
  AffineExpr(
    BigInt(c) * a.constant, Dict{Int,BigInt}(k => BigInt(c) * v for (k, v) in a.coeffs)
  )
end

Base.:*(a::AffineExpr, c::Integer) = c * a

function Base.:(==)(a::AffineExpr, b::AffineExpr)
  a.constant == b.constant && a.coeffs == b.coeffs
end
Base.:(==)(a::AffineExpr, c::Integer) = is_determined(a) && a.constant == c
Base.:(==)(c::Integer, a::AffineExpr) = a == c

Base.:+(a::AffineExpr, c::Integer) = a + AffineExpr(c)
Base.:+(c::Integer, a::AffineExpr) = a + c
Base.:-(a::AffineExpr, c::Integer) = a - AffineExpr(c)
Base.:-(c::Integer, a::AffineExpr) = AffineExpr(c) - a

"""Alternating sum ``e_1 - e_2 + e_3 - \\cdots`` of the entries."""
function _alternating_sum(entries::AbstractVector{AffineExpr})
  total = AffineExpr(0)
  for (i, entry) in enumerate(entries)
    total = isodd(i) ? total + entry : total - entry
  end
  total
end

function _alternating_sum(M::AbstractMatrix{AffineExpr}, row::Int, dim::Int)
  _alternating_sum(@view M[row, 1:(dim + 1)])
end

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, e::AffineExpr)
  is_determined(e) && return print(io, e.constant)
  leading = e.constant == 0
  leading || print(io, e.constant)
  for (var_id, coeff) in sort(collect(e.coeffs); by=first)
    if leading
      # First printed term: sign is attached, "1" is suppressed.
      coeff == 1 || print(io, coeff == -1 ? "-" : "$coeff * ")
      leading = false
    else
      print(io, coeff > 0 ? " + " : " - ")
      abs(coeff) == 1 || print(io, abs(coeff), " * ")
    end
    print(io, "x_$var_id")
  end
end

function Base.show(io::IO, H::Cohomology{AffineExpr})
  parts = [
    "H$(_superscript(i)) = $(sprint(show, H[i]))" for
    i in 0:(H.max_degree) if !is_zero_expr(H[i])
  ]
  print(io, isempty(parts) ? "H* = 0" : join(parts, "\n"))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Linear constraint propagation
# ═══════════════════════════════════════════════════════════════════════════════

"""
Substitute ``x_{\\mathrm{var\\_id}} = \\mathrm{replacement}`` in every entry
of the `AffineExpr` array ``M``, in place.  Returns `true` when at least one
entry changed.
"""
function _substitute_var!(
  M::AbstractArray{AffineExpr}, var_id::Int, replacement::AffineExpr
)
  changed = false
  for i in eachindex(M)
    entry = M[i]
    haskey(entry.coeffs, var_id) || continue
    factor = entry.coeffs[var_id]

    # entry = (rest) + factor * x_{var_id}  ⟶  (rest) + factor * replacement
    new_constant = entry.constant + factor * replacement.constant
    new_coeffs = copy(entry.coeffs)
    delete!(new_coeffs, var_id)
    for (var, coeff) in replacement.coeffs
      new_coeffs[var] = get(new_coeffs, var, BigInt(0)) + factor * coeff
      new_coeffs[var] == 0 && delete!(new_coeffs, var)
    end
    M[i] = AffineExpr(new_constant, new_coeffs)
    changed = true
  end
  changed
end

"""
Solve the linear equation ``\\mathrm{expr} = 0`` for one of its variables:
the candidate with the smallest id whose coefficient divides the constant
and all other coefficients, so that the substitution stays integral.
Returns `var_id => solution`, or `nothing` when no variable qualifies.

With `among` given, only variables from that set are candidates (the
solution itself still involves all variables of the equation).

The default choice depends on the equation only, never on the array the
solution is applied to: fixed-point loops that apply one equation to
several arrays must eliminate the *same* variable everywhere, otherwise two
arrays can trade variables back and forth without ever converging.
"""
function _solve_for_variable(expr::AffineExpr; among::Union{Nothing,Set{Int}}=nothing)
  candidates = if among === nothing
    collect(keys(expr.coeffs))
  else
    [var for var in keys(expr.coeffs) if var in among]
  end

  for var_id in sort!(candidates)
    coeff = expr.coeffs[var_id]
    rest_coeffs = copy(expr.coeffs)
    delete!(rest_coeffs, var_id)

    # Integrality check: x_{var_id} = -(constant + rest) / coeff must have
    # integer coefficients, since every variable is an integer.
    expr.constant % coeff == 0 || continue
    all(v % coeff == 0 for (_, v) in rest_coeffs) || continue

    solution_constant = -(expr.constant ÷ coeff)
    solution_coeffs = Dict{Int,BigInt}(k => -(v ÷ coeff) for (k, v) in rest_coeffs)
    filter!(p -> p.second != 0, solution_coeffs)
    return var_id => AffineExpr(solution_constant, solution_coeffs)
  end
  nothing
end

"""
Given a linear equation ``\\mathrm{expr} = 0`` in the symbolic variables,
eliminate one variable (chosen by `_solve_for_variable`, i.e.
consistently across arrays) from the `AffineExpr` array ``M``, in place.
Returns `true` when an entry of ``M`` changed.
"""
function _apply_equation!(M::AbstractArray{AffineExpr}, expr::AffineExpr)
  solution = _solve_for_variable(expr)
  solution === nothing && return false
  _substitute_var!(M, solution.first, solution.second)
end

"""
Like [`_apply_equation!`](@ref), but only eliminates a variable that occurs
in ``M`` itself; variables of the equation that do not occur in ``M`` are
kept as parameters.  This is what the entry-based LES solvers need: they
must solve for their own fresh output variables, not for the variables of
the input terms.
"""
function _apply_equation_in_vars!(M::AbstractArray{AffineExpr}, expr::AffineExpr)
  isempty(expr.coeffs) && return false

  array_vars = Set{Int}()
  for entry in M, var in keys(entry.coeffs)
    push!(array_vars, var)
  end

  solution = _solve_for_variable(expr; among=array_vars)
  solution === nothing && return false
  _substitute_var!(M, solution.first, solution.second)
end

"""
Renumber the symbolic variables of an `AffineExpr` array to contiguous ids
``x_0, x_1, \\ldots`` so the output does not reflect the internal counter.
"""
function _renumber_variables!(M::AbstractArray{AffineExpr})
  old_ids = Set{Int}()
  for entry in M, var in keys(entry.coeffs)
    push!(old_ids, var)
  end
  isempty(old_ids) && return M

  mapping = Dict{Int,Int}(
    old_id => new_id - 1 for (new_id, old_id) in enumerate(sort!(collect(old_ids)))
  )
  for i in eachindex(M)
    entry = M[i]
    isempty(entry.coeffs) && continue
    M[i] = AffineExpr(
      entry.constant, Dict{Int,BigInt}(mapping[k] => v for (k, v) in entry.coeffs)
    )
  end
  M
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Interval propagation
# ═══════════════════════════════════════════════════════════════════════════════

_infeasible_intervals() = error(
  "inconsistent long exact sequence constraints: the input data admits no " *
  "nonnegative solution (is the section regular, with smooth zero locus of " *
  "the expected dimension?)",
)

"""
    _propagate_intervals!(system::AbstractArray{AffineExpr}) -> Bool

Interval (bound-consistency) propagation on a system of nonnegative
quantities.

Every element of `system` represents a nonnegative integer (a cohomology
dimension or a connecting-map rank), and every symbolic variable is itself a
nonnegative integer.  Each element ``c + \\sum_j k_j x_j \\ge 0`` bounds each of
its variables in terms of the current intervals of the others; the bounds are
tightened to a fixed point, and a variable whose interval collapses to a
point is substituted throughout `system`.  Raises when an interval empties,
which means the input data was inconsistent.

Returns `true` when at least one variable was pinned.  This recovers the
information the equation-only solvers discard: for example the entries
``(x - 20, 21 - x, 20 - x + y)`` pin ``x = 20`` and then bound ``y``.
"""
function _propagate_intervals!(system::AbstractArray{AffineExpr})
  # Fast path: fully determined systems (the common case) need no sweeps,
  # only the nonnegativity sanity check.
  if all(is_determined, system)
    all(constraint -> constraint.constant >= 0, system) || _infeasible_intervals()
    return false
  end

  lower = Dict{Int,BigInt}()  # absent: 0 (every variable is a nonnegative integer)
  upper = Dict{Int,BigInt}()  # absent: unbounded

  pinned_any = false
  # The sweep cap is a backstop: bounds tighten monotonically, so feasible
  # systems reach their fixed point long before it (stopping early is sound).
  for _ in 1:64
    changed = false

    for constraint in system
      if isempty(constraint.coeffs)
        constraint.constant >= 0 || _infeasible_intervals()
        continue
      end
      for (var, coeff) in constraint.coeffs
        # From c + Σⱼ kⱼ xⱼ ≥ 0, isolate the term of `var`:
        #     coeff * x_var ≥ -(c + Σ_{j≠var} kⱼ xⱼ) ≥ -residual,
        # where `residual` is the box maximum of c + Σ_{j≠var} kⱼ xⱼ
        # (upper bounds for positive kⱼ, lower bounds for negative kⱼ).
        # Skip when that maximum is unbounded.
        residual = constraint.constant
        bounded = true
        for (other, other_coeff) in constraint.coeffs
          other == var && continue
          if other_coeff > 0
            other_upper = get(upper, other, nothing)
            other_upper === nothing && (bounded=false; break)
            residual += other_coeff * other_upper
          else
            residual += other_coeff * get(lower, other, BigInt(0))
          end
        end
        bounded || continue

        if coeff < 0
          # (-coeff) * x_var ≤ residual, so x_var ≤ ⌊residual / (-coeff)⌋.
          new_upper = fld(residual, -coeff)
          if new_upper < get(upper, var, new_upper + 1)
            upper[var] = new_upper
            changed = true
            new_upper < get(lower, var, BigInt(0)) && _infeasible_intervals()
          end
        else
          # coeff * x_var ≥ -residual, so x_var ≥ ⌈-residual / coeff⌉.
          new_lower = cld(-residual, coeff)
          if new_lower > get(lower, var, BigInt(0))
            lower[var] = new_lower
            changed = true
            new_lower > get(upper, var, new_lower) && _infeasible_intervals()
          end
        end
      end
    end

    # A collapsed interval turns the variable into a known constant.
    for var in collect(keys(upper))
      if get(lower, var, BigInt(0)) == upper[var]
        _substitute_var!(system, var, AffineExpr(upper[var]))
        delete!(upper, var)
        delete!(lower, var)
        changed = true
        pinned_any = true
      end
    end

    changed || break
  end
  pinned_any
end

"""
Impose the vanishing ``\\mathrm{H}^k = 0`` for ``k > d`` on a symbolic cohomology
vector (indexed by degree ``0, \\ldots, \\mathrm{length} - 1``), eliminating
symbolic variables where possible, and truncate to degrees ``0, \\ldots, d``.
"""
function _truncate_cohomology!(entries::Vector{AffineExpr}, d::Int)
  for k in (d + 2):length(entries)
    is_zero_expr(entries[k]) || _apply_equation!(entries, entries[k])
  end
  entries[1:(d + 1)]
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Short exact sequence solver (BigInt-valued)
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compute bounds on the connecting-map ranks of a short exact sequence
``0 \\to \\mathcal{A} \\to \\mathcal{B} \\to \\mathcal{C} \\to 0`` given integer cohomology dimensions for
``\\mathcal{A}`` and ``\\mathcal{B}``.

Writing ``δ_i = \\mathrm{rank}(\\mathrm{H}^i(X, \\mathcal{C}) \\to \\mathrm{H}^{i+1}(X, \\mathcal{A}))``, exactness gives
``c_i = b_i - a_i + δ_{i-1} + δ_i`` with ``0 \\le δ_i \\le a_{i+1}`` and
``δ_{i-1} + δ_i \\ge a_i - b_i``.  The bounds are propagated
forward and backward until convergence.

Returns `(lower, upper)` bound vectors of length `d + 2` indexed for ``δ_{-1}, …, δ_d``.
"""
function _ses_delta_bounds(a_vals::Vector{BigInt}, b_vals::Vector{BigInt}, d::Int)
  # Every rank has a finite bound from the start, so no infinity sentinel is
  # needed: δ_{-1} = δ_d = 0 (there is no H^{-1}(C) and no H^{d+1}(A)), and
  # δ_i ≤ a_{i+1} since the connecting map lands in H^{i+1}(A).
  upper = [BigInt(0); a_vals[2:end]; BigInt(0)]

  # δ_i ≥ a_{i+1} - b_{i+1}: exactness at H^{i+1}(A) makes the kernel of
  # H^{i+1}(A) → H^{i+1}(B) the image of the connecting map.
  lower = zeros(BigInt, d + 2)
  for i in 0:(d - 1)
    lower[i + 2] = max(lower[i + 2], a_vals[i + 2] - b_vals[i + 2])
  end

  # Each pass can only push information one slot along the sequence, so
  # d + 2 passes reach the fixed point.
  for _ in 1:(d + 2)
    changed = false

    # c_i ≥ 0 gives δ_{i-1} + δ_i ≥ a_i - b_i: a small upper bound on either
    # rank forces a lower bound on the other.
    for i in 0:d
      needed = a_vals[i + 1] - b_vals[i + 1]
      for (target, partner) in ((i + 2, i + 1), (i + 1, i + 2))
        new_lb = needed - upper[partner]
        if new_lb > lower[target]
          lower[target] = new_lb
          changed = true
        end
      end
    end

    # Clamp: an interval that would empty is truncated (the caller treats
    # lower == upper as determined and anything else as open).
    for j in 1:(d + 2)
      lower[j] > upper[j] && (lower[j] = upper[j])
    end

    !changed && break
  end

  (lower, upper)
end

"""
    solve_ses_cohomology(a::Cohomology{BigInt}, b::Cohomology{BigInt})
      -> (Cohomology{BigInt}, Bool)

Given ``\\mathrm{H}^\\bullet(X, \\mathcal{A})`` and ``\\mathrm{H}^\\bullet(X, \\mathcal{B})`` from a short exact sequence
``0 \\to \\mathcal{A} \\to \\mathcal{B} \\to \\mathcal{C} \\to 0``, determine ``\\mathrm{H}^\\bullet(X, \\mathcal{C})`` via the long
exact sequence.

Returns the cohomology of ``\\mathcal{C}`` together with a Boolean flag
`determined`, which is `true` when the long exact sequence uniquely
determines all cohomology groups.  When the flag is `false` the entries are
lower bounds and must not be used as values; use the symbolic solvers
instead.
"""
function solve_ses_cohomology(a::Cohomology{BigInt}, b::Cohomology{BigInt})
  d = a.max_degree
  b.max_degree == d || throw(
    ArgumentError(
      "solve_ses_cohomology: a and b must have equal max_degree, got $d and $(b.max_degree)"
    ),
  )

  a_vals = BigInt[a[i] for i in 0:d]
  b_vals = BigInt[b[i] for i in 0:d]
  (lower, upper) = _ses_delta_bounds(a_vals, b_vals, d)

  # c_i = b_i - a_i + δ_{i-1} + δ_i, using the lower bounds for the δ's.
  entries = BigInt[b_vals[i + 1] - a_vals[i + 1] + lower[i + 1] + lower[i + 2] for i in 0:d]
  (Cohomology{BigInt}(entries, d), lower == upper)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul filtration solver (BigInt-valued)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_koszul_filtration(
      koszul_cohos::Vector{Cohomology{BigInt}},
      dim_zero_locus::Int
    ) -> (Cohomology{BigInt}, Bool)

Given ``\\mathrm{H}^\\bullet(X, \\mathcal{K}_i)`` for the Koszul terms
``\\mathcal{K}_i = \\mathcal{F} \\otimes \\wedge^i \\mathcal{E}^\\vee`` (``i = 0, 1, \\ldots, r``),
compute ``\\mathrm{H}^\\bullet(Z, \\mathcal{F}|_Z)`` by iterating through the short exact sequences
arising from the Koszul filtration:

```math
\\mathcal{C}_r = \\mathcal{K}_r, \\qquad
0 \\to \\mathcal{C}_{j+1} \\to \\mathcal{K}_j \\to \\mathcal{C}_j \\to 0
\\qquad (j = r-1, \\ldots, 0)
```

where ``\\mathcal{C}_0 \\cong \\mathcal{F}|_Z`` (shifted to the zero locus dimension).

Returns the cohomology of ``\\mathcal{F}|_Z`` together with a `determined`
flag.  When the chain of short exact sequences
leaves some connecting maps open, a symbolic re-solve with the vanishing
``\\mathrm{H}^k(Z, \\mathcal{F}|_Z) = 0`` for ``k > \\dim Z`` and the exact Euler characteristic
is attempted before giving up.
"""
function solve_koszul_filtration(
  koszul_cohos::Vector{Cohomology{BigInt}},
  dim_zero_locus::Int,
)
  r = length(koszul_cohos) - 1  # rank of the defining bundle

  # Chain from the top of the filtration downwards: C_r = K_r, and each step
  # solves 0 → C_{j+1} → K_j → C_j → 0 for its cokernel.
  current = koszul_cohos[r + 1]
  all_determined = true
  for j in (r - 1):-1:0
    (current, determined) = solve_ses_cohomology(current, koszul_cohos[j + 1])
    all_determined &= determined
  end

  result = Cohomology{BigInt}(BigInt[current[i] for i in 0:dim_zero_locus], dim_zero_locus)
  all_determined && return (result, true)

  # Some intermediate SES was undetermined, but the final output may still be
  # pinned down: re-solve symbolically over the full ambient range, then apply
  # the vanishing H^k(Z, F|_Z) = 0 for k > dim Z and the exact Euler
  # characteristic χ(F|_Z) = Σ_i (-1)^i χ(K_i).
  var_counter = Ref(0)
  dim_ambient = koszul_cohos[1].max_degree
  sym = solve_koszul_filtration_symbolic(koszul_cohos, dim_ambient, var_counter)
  entries = _truncate_cohomology!(AffineExpr[sym[k] for k in 0:dim_ambient], dim_zero_locus)

  chi_exact = _alternating_euler_characteristic(koszul_cohos)
  _apply_equation!(entries, _alternating_sum(entries) - AffineExpr(chi_exact))

  if all(is_determined, entries)
    return (Cohomology{BigInt}(_determined_bigints(entries), dim_zero_locus), true)
  end

  (result, false)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Short exact sequence solvers (symbolic, δ-based)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_ses_cohomology_symbolic(a, b, var_counter) -> Cohomology{AffineExpr}

Symbolic version of [`solve_ses_cohomology`](@ref): solve
``0 \\to \\mathcal{A} \\to \\mathcal{B} \\to \\mathcal{C} \\to 0`` for ``\\mathrm{H}^\\bullet(X, \\mathcal{C})``, introducing a fresh
symbolic variable for every connecting-map rank that is not forced.

For `Cohomology{BigInt}` inputs the connecting-map bounds are
propagated first and variables are only introduced where a gap remains.
For symbolic inputs each rank ``δ_i`` is zero when ``\\mathrm{H}^{i+1}(X, \\mathcal{A}) = 0``,
equal to ``a_{i+1}`` when ``\\mathrm{H}^{i+1}(X, \\mathcal{B}) = 0`` (the connecting map is then
surjective), and a fresh variable otherwise.

The `var_counter` is advanced for each new variable; pass a shared counter
across multiple calls to keep variable ids globally distinct.
"""
function solve_ses_cohomology_symbolic(
  a::Cohomology{BigInt}, b::Cohomology{BigInt},
  var_counter::Ref{Int},
)
  d = a.max_degree
  b.max_degree == d || throw(
    ArgumentError(
      "solve_ses_cohomology_symbolic: a and b must have equal max_degree, got $d and $(b.max_degree)"
    ),
  )

  a_vals = BigInt[a[i] for i in 0:d]
  b_vals = BigInt[b[i] for i in 0:d]
  (lower, upper) = _ses_delta_bounds(a_vals, b_vals, d)

  # δ_j = lb_j exactly when the bounds meet, and lb_j plus a fresh
  # nonnegative variable otherwise.
  δ = AffineExpr[
    if lower[j] == upper[j]
      AffineExpr(lower[j])
    else
      AffineExpr(lower[j]) + _fresh_variable(var_counter)
    end
    for j in 1:(d + 2)
  ]

  entries = AffineExpr[
    AffineExpr(b_vals[i + 1] - a_vals[i + 1]) + δ[i + 1] + δ[i + 2] for i in 0:d
  ]
  Cohomology{AffineExpr}(entries, d)
end

function solve_ses_cohomology_symbolic(
  a::Cohomology{AffineExpr}, b::Cohomology{AffineExpr},
  var_counter::Ref{Int},
)
  d = a.max_degree
  b.max_degree == d || throw(
    ArgumentError(
      "solve_ses_cohomology_symbolic: a and b must have equal max_degree, got $d and $(b.max_degree)"
    ),
  )

  # Fully determined inputs: BigInt bound propagation is sharper.
  if all(is_determined(a[i]) for i in 0:d) && all(is_determined(b[i]) for i in 0:d)
    return solve_ses_cohomology_symbolic(
      Cohomology{BigInt}(BigInt[a[i].constant for i in 0:d], d),
      Cohomology{BigInt}(BigInt[b[i].constant for i in 0:d], d),
      var_counter,
    )
  end

  δ = Vector{AffineExpr}(undef, d + 2)
  δ[1] = AffineExpr(0)      # δ_{-1} = 0
  δ[d + 2] = AffineExpr(0)  # δ_d = 0
  for i in 0:(d - 1)
    δ[i + 2] = if is_zero_expr(a[i + 1])
      AffineExpr(0)         # nothing to land in
    elseif is_zero_expr(b[i + 1])
      a[i + 1]              # surjective connecting map
    else
      _fresh_variable(var_counter)
    end
  end

  entries = AffineExpr[b[i] - a[i] + δ[i + 1] + δ[i + 2] for i in 0:d]
  Cohomology{AffineExpr}(entries, d)
end

function solve_ses_cohomology_symbolic(
  a::Cohomology{AffineExpr}, b::Cohomology{BigInt},
  var_counter::Ref{Int},
)
  d = b.max_degree
  b_affine = Cohomology{AffineExpr}(_as_affine(BigInt[b[i] for i in 0:d]), d)
  solve_ses_cohomology_symbolic(a, b_affine, var_counter)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul filtration solver (symbolic, δ-based)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_koszul_filtration_symbolic(
      koszul_cohos::Vector{Cohomology{BigInt}},
      dim_zero_locus::Int,
      var_counter::Ref{Int},
    ) -> Cohomology{AffineExpr}

Symbolic version of [`solve_koszul_filtration`](@ref): introduces fresh
symbolic variables for each undetermined connecting-map rank rather than
using lower bounds.
"""
function solve_koszul_filtration_symbolic(
  koszul_cohos::Vector{Cohomology{BigInt}},
  dim_zero_locus::Int,
  var_counter::Ref{Int},
)
  r = length(koszul_cohos) - 1

  if r == 0
    entries = AffineExpr[AffineExpr(koszul_cohos[1][i]) for i in 0:dim_zero_locus]
    return Cohomology{AffineExpr}(entries, dim_zero_locus)
  end

  # C_r = K_r, then 0 → C_{j+1} → K_j → C_j → 0 for j = r-1, …, 0.  The first
  # step consumes K_r and K_{r-1} at once (both BigInt-valued), which is why the
  # loop starts at j = r - 2.
  current = solve_ses_cohomology_symbolic(
    koszul_cohos[r + 1], koszul_cohos[r], var_counter
  )
  for j in (r - 2):-1:0
    current = solve_ses_cohomology_symbolic(current, koszul_cohos[j + 1], var_counter)
  end

  entries = AffineExpr[current[i] for i in 0:dim_zero_locus]
  Cohomology{AffineExpr}(entries, dim_zero_locus)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Entry-based LES solvers
#
#  Alternative to the δ-based symbolic solvers above: introduce a fresh
#  symbolic variable for each *unknown* entry of the long exact sequence and
#  derive linear equations from the alternating-sum condition on the maximal
#  nonzero segments.  This matches the approach in Macaulay2's
#  `shortExactSequenceCoker`.  The partial alternating sums additionally
#  provide inequalities, consumed by `_propagate_intervals!`.
# ═══════════════════════════════════════════════════════════════════════════════

"""
Alternating-sum equations extracted from a long exact sequence.

The sequence is split at entries that are identically zero; exactness forces
the alternating sum of every maximal nonzero segment to vanish.  An entry
that is symbolic but happens to take the value zero only merges two segments,
and the sum of two valid equations is still valid, so the extracted equations
hold in all cases.
"""
function _les_equations(les::Vector{AffineExpr})
  equations = AffineExpr[]
  segment = AffineExpr[]
  for entry in les
    if is_zero_expr(entry)
      isempty(segment) || push!(equations, _alternating_sum(segment))
      segment = AffineExpr[]
    else
      push!(segment, entry)
    end
  end
  isempty(segment) || push!(equations, _alternating_sum(segment))
  equations
end

"""
Partial alternating-sum inequalities of a long exact sequence with zero
flanks: writing ``r_k`` for the rank of the map out of slot ``k``, exactness
gives ``s_k = r_{k-1} + r_k``, so the partial sums satisfy

```math
(-1)^{k-1} \\sum_{i \\le k} (-1)^{i-1} s_i = r_k \\ge 0 .
```

Returns one expression per proper partial sum (the full sum is the usual
equality and is omitted), each of which is a nonnegative quantity.  These
encode every linear consequence of exactness, in particular the bounds
``r_k \\le \\min(s_k, s_{k+1})`` that the equation-only solvers discard.
Unlike `_les_equations` this needs no segmentation: the whole interleaved
sequence is exact with zero flanks.
"""
function _les_inequalities(les::Vector{AffineExpr})
  inequalities = AffineExpr[]
  partial_sum = AffineExpr(0)
  for k in 1:(length(les) - 1)
    partial_sum = isodd(k) ? partial_sum + les[k] : partial_sum - les[k]
    rank_k = isodd(k) ? partial_sum : -partial_sum
    is_determined(rank_k) || push!(inequalities, rank_k)
  end
  inequalities
end

"""
Interleave the cohomologies of a short exact sequence
``0 \\to \\mathcal{A} \\to \\mathcal{B} \\to \\mathcal{C} \\to 0`` into its long exact sequence
``(a_0, b_0, c_0, a_1, b_1, c_1, \\ldots)``.
"""
function _les_interleave(
  a::Vector{AffineExpr}, b::Vector{AffineExpr}, c::Vector{AffineExpr}
)
  n = length(a)
  les = Vector{AffineExpr}(undef, 3n)
  for i in 1:n
    les[3i - 2] = a[i]
    les[3i - 1] = b[i]
    les[3i] = c[i]
  end
  les
end

function _determined_les_cokernel(a::Vector{BigInt}, b::Vector{BigInt})
  d = length(a) - 1
  c, determined = solve_ses_cohomology(Cohomology{BigInt}(a, d), Cohomology{BigInt}(b, d))
  determined || return nothing
  AffineExpr[AffineExpr(c[i]) for i in 0:d]
end

"""
    les_cokernel(a, b, var_counter; inequalities=nothing) -> Vector{AffineExpr}

Given ``\\mathrm{H}^\\bullet(X, \\mathcal{A})`` and ``\\mathrm{H}^\\bullet(X, \\mathcal{B})`` from a short exact sequence
``0 \\to \\mathcal{A} \\to \\mathcal{B} \\to \\mathcal{C} \\to 0``, compute ``\\mathrm{H}^\\bullet(X, \\mathcal{C})``.

Creates a fresh symbolic variable for each ``\\mathrm{H}^i(X, \\mathcal{C})``, then eliminates as
many as possible using the alternating-sum equations of the long exact
sequence (see `_les_equations`). Fully determined input is
delegated to the bound-propagation solver first.

When an `inequalities` vector is passed, the partial alternating-sum
inequalities of the sequence (see `_les_inequalities`) are appended to it,
in terms of the surviving variables, for later interval propagation.
"""
function les_cokernel(
  a::Vector{AffineExpr}, b::Vector{AffineExpr},
  var_counter::Ref{Int};
  inequalities::Union{Nothing,Vector{AffineExpr}}=nothing,
)
  n = length(a)
  length(b) == n || throw(
    ArgumentError("les_cokernel: a and b must have equal length, got $n and $(length(b))")
  )

  # H^*(A) = 0 makes B → C an isomorphism on cohomology.
  all(is_zero_expr, a) && return copy(b)

  # Fully determined input: bound propagation is sharper than the segment
  # equations, use it whenever it determines the answer.
  if all(is_determined, a) && all(is_determined, b)
    determined_result = _determined_les_cokernel(
      _determined_bigints(a), _determined_bigints(b)
    )
    determined_result !== nothing && return determined_result
  end

  # One fresh variable per unknown entry of C, then eliminate through the
  # alternating-sum equations of the interleaved sequence.
  c = AffineExpr[_fresh_variable(var_counter) for _ in 1:n]
  for eq in _les_equations(_les_interleave(a, b, c))
    _apply_equation_in_vars!(c, eq)
  end
  if inequalities !== nothing
    # Harvest after elimination so the inequalities are stated in the
    # variables that actually survive.
    append!(inequalities, _les_inequalities(_les_interleave(a, b, c)))
  end
  c
end

function les_cokernel(
  a::Vector{AffineExpr}, b::Vector{BigInt}, var_counter::Ref{Int};
  inequalities::Union{Nothing,Vector{AffineExpr}}=nothing,
)
  les_cokernel(a, _as_affine(b), var_counter; inequalities)
end

"""
    les_kernel(b, c, var_counter; inequalities=nothing) -> Vector{AffineExpr}

Given ``\\mathrm{H}^\\bullet(X, \\mathcal{B})`` and ``\\mathrm{H}^\\bullet(X, \\mathcal{C})`` from a short exact sequence
``0 \\to \\mathcal{A} \\to \\mathcal{B} \\to \\mathcal{C} \\to 0``, compute ``\\mathrm{H}^\\bullet(X, \\mathcal{A})``.

Dual to [`les_cokernel`](@ref): the unknowns sit in the ``\\mathcal{A}``-slots of the
long exact sequence.
"""
function les_kernel(
  b::Vector{AffineExpr}, c::Vector{AffineExpr},
  var_counter::Ref{Int};
  inequalities::Union{Nothing,Vector{AffineExpr}}=nothing,
)
  n = length(b)
  length(c) == n || throw(
    ArgumentError("les_kernel: b and c must have equal length, got $n and $(length(c))")
  )

  # H^*(C) = 0 makes A → B an isomorphism on cohomology.
  all(is_zero_expr, c) && return copy(b)

  a = AffineExpr[_fresh_variable(var_counter) for _ in 1:n]
  for eq in _les_equations(_les_interleave(a, b, c))
    _apply_equation_in_vars!(a, eq)
  end
  if inequalities !== nothing
    append!(inequalities, _les_inequalities(_les_interleave(a, b, c)))
  end
  a
end

"""
    long_exact_sequence_cokernel(terms, var_counter; inequalities=nothing)
      -> Vector{AffineExpr}

Given cohomology of terms `[K_r, K_{r-1}, …, K_0]` (reversed Koszul order),
iteratively apply [`les_cokernel`](@ref) to compute the final cokernel:

```math
\\mathcal{C}_r = \\mathcal{K}_r, \\qquad
0 \\to \\mathcal{C}_{j+1} \\to \\mathcal{K}_j \\to \\mathcal{C}_j \\to 0
\\qquad (j = r-1, \\ldots, 0)
```
"""
function long_exact_sequence_cokernel(
  terms::Vector{Vector{T}}, var_counter::Ref{Int};
  inequalities::Union{Nothing,Vector{AffineExpr}}=nothing,
) where {T<:Union{BigInt,AffineExpr}}
  current = _as_affine(terms[1])
  for j in 2:length(terms)
    current = les_cokernel(current, terms[j], var_counter; inequalities)
  end
  current
end
