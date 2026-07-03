# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul.jl — Long exact sequence algebra for Koszul resolutions
#
#  Pure algebra module: solves long exact sequences arising from
#  Koszul-type filtrations.  Knows nothing about flag varieties or bundles.
#
#  Given the dimension-valued cohomology of two terms of a short exact
#  sequence, computes the cohomology of the third term as far as exactness
#  determines it.  Undetermined connecting-map ranks are either reported
#  (numeric solvers, which return a `determined` flag) or turned into
#  symbolic variables (the `AffineExpr`-valued solvers).
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

"""Create a fresh symbolic variable and advance the counter."""
function _fresh_variable(var_counter::Ref{Int})
  x = symbolic_variable(var_counter[])
  var_counter[] += 1
  x
end

"""Check whether the expression is fully determined (no symbolic variables)."""
is_determined(e::AffineExpr) = isempty(e.coeffs)

"""Check whether the expression is identically zero."""
is_zero_expr(e::AffineExpr) = e.constant == 0 && isempty(e.coeffs)

_determined_bigints(entries::Vector{AffineExpr}) = BigInt[e.constant for e in entries]

_as_affine(v::Vector{AffineExpr}) = copy(v)
_as_affine(v::Vector{BigInt}) = AffineExpr.(v)

# ─── Arithmetic ──────────────────────────────────────────────────────────────

function _merge_coeffs(op, a::Dict{Int,BigInt}, b::Dict{Int,BigInt})
  result = copy(a)
  for (k, v) in b
    result[k] = op(get(result, k, BigInt(0)), v)
    result[k] == 0 && delete!(result, k)
  end
  result
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

function _alternating_sum(entries::AbstractVector{AffineExpr})
  total = AffineExpr(0)
  add_term = true
  for entry in entries
    total = add_term ? total + entry : total - entry
    add_term = !add_term
  end
  total
end

function _alternating_sum(mat::AbstractMatrix{AffineExpr}, row::Int, dim::Int)
  _alternating_sum(@view mat[row, 1:(dim + 1)])
end

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, e::AffineExpr)
  if isempty(e.coeffs)
    print(io, e.constant)
    return nothing
  end
  sorted_vars = sort(collect(e.coeffs); by=first)
  first_term = true
  if e.constant != 0
    print(io, e.constant)
    first_term = false
  end
  for (var_id, coeff) in sorted_vars
    if first_term
      if coeff == 1
        print(io, "x_$var_id")
      elseif coeff == -1
        print(io, "-x_$var_id")
      else
        print(io, "$coeff * x_$var_id")
      end
      first_term = false
    else
      if coeff > 0
        coeff == 1 ? print(io, " + x_$var_id") : print(io, " + $coeff * x_$var_id")
      else
        coeff == -1 ? print(io, " - x_$var_id") : print(io, " - $(-coeff) * x_$var_id")
      end
    end
  end
end

function Base.show(io::IO, H::Cohomology{AffineExpr})
  parts = String[]
  for i in 0:(H.dim_variety)
    v = H[i]
    is_zero_expr(v) && continue
    push!(parts, "H$(_superscript(i)) = $(sprint(show, v))")
  end
  if isempty(parts)
    print(io, "H* = 0")
  else
    print(io, join(parts, "\n"))
  end
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
    e = M[i]
    haskey(e.coeffs, var_id) || continue
    c = e.coeffs[var_id]
    new_constant = e.constant + c * replacement.constant
    new_coeffs = copy(e.coeffs)
    delete!(new_coeffs, var_id)
    for (k, v) in replacement.coeffs
      new_coeffs[k] = get(new_coeffs, k, BigInt(0)) + c * v
      new_coeffs[k] == 0 && delete!(new_coeffs, k)
    end
    M[i] = AffineExpr(new_constant, new_coeffs)
    changed = true
  end
  changed
end

"""
Solve the linear equation ``\\mathrm{expr} = 0`` for one of its variables:
the one with the smallest id whose coefficient divides all other
coefficients (so the substitution stays integral).  Returns
`var_id => solution`, or `nothing` when no variable qualifies.

The choice depends on the equation only, never on the array the solution is
applied to: fixed-point loops that apply one equation to several arrays
must eliminate the *same* variable everywhere, otherwise two arrays can
trade variables back and forth without ever converging.
"""
function _solve_for_variable(expr::AffineExpr)
  for var_id in sort!(collect(keys(expr.coeffs)))
    coeff = expr.coeffs[var_id]
    rest_coeffs = copy(expr.coeffs)
    delete!(rest_coeffs, var_id)
    expr.constant % coeff == 0 || continue
    all(v % coeff == 0 for (_, v) in rest_coeffs) || continue

    sub_const = -(expr.constant ÷ coeff)
    sub_coeffs = Dict{Int,BigInt}(k => -(v ÷ coeff) for (k, v) in rest_coeffs)
    filter!(p -> p.second != 0, sub_coeffs)
    return var_id => AffineExpr(sub_const, sub_coeffs)
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
  for e in M, v in keys(e.coeffs)
    push!(array_vars, v)
  end

  restricted_coeffs = Dict{Int,BigInt}(
    v => c for (v, c) in expr.coeffs if v in array_vars
  )
  isempty(restricted_coeffs) && return false

  # Solve, but only allow eliminating a variable present in M: keep the
  # full equation (all variables) as the solution content.
  for var_id in sort!(collect(keys(restricted_coeffs)))
    coeff = expr.coeffs[var_id]
    rest_coeffs = copy(expr.coeffs)
    delete!(rest_coeffs, var_id)
    expr.constant % coeff == 0 || continue
    all(v % coeff == 0 for (_, v) in rest_coeffs) || continue

    sub_const = -(expr.constant ÷ coeff)
    sub_coeffs = Dict{Int,BigInt}(k => -(v ÷ coeff) for (k, v) in rest_coeffs)
    filter!(p -> p.second != 0, sub_coeffs)
    return _substitute_var!(M, var_id, AffineExpr(sub_const, sub_coeffs))
  end
  false
end

"""
Renumber the symbolic variables of an `AffineExpr` array to contiguous ids
``x_0, x_1, \\ldots`` so the output does not reflect the internal counter.
"""
function _renumber_variables!(M::AbstractArray{AffineExpr})
  old_ids = Set{Int}()
  for e in M, v in keys(e.coeffs)
    push!(old_ids, v)
  end
  isempty(old_ids) && return M

  mapping = Dict{Int,Int}(
    old_id => new_id - 1 for (new_id, old_id) in enumerate(sort!(collect(old_ids)))
  )
  for i in eachindex(M)
    e = M[i]
    isempty(e.coeffs) && continue
    M[i] = AffineExpr(
      e.constant, Dict{Int,BigInt}(mapping[k] => v for (k, v) in e.coeffs)
    )
  end
  M
end

"""
Impose the vanishing ``H^k = 0`` for ``k > d`` on a symbolic cohomology
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
#  Short exact sequence solver (numeric)
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compute bounds on the connecting-map ranks of a short exact sequence
``0 \\to A \\to B \\to C \\to 0`` given numeric cohomology values for
``A`` and ``B``.

Writing ``δ_i = \\mathrm{rank}(H^i(C) \\to H^{i+1}(A))``, exactness gives
``c_i = b_i - a_i + δ_{i-1} + δ_i`` with ``0 \\le δ_i \\le a_{i+1}`` and
``δ_{i-1} + δ_i \\ge a_i - b_i``.  The bounds are propagated
forward and backward until convergence.

Returns `(lb, ub)`, vectors of length `d + 2` indexed for ``δ_{-1}, …, δ_d``.
"""
function _ses_delta_bounds(a_vals::Vector{BigInt}, b_vals::Vector{BigInt}, d::Int)
  lb = zeros(BigInt, d + 2)
  ub = fill(BigInt(10)^18, d + 2)

  ub[1] = BigInt(0)      # δ_{-1} = 0
  ub[d + 2] = BigInt(0)  # δ_d = 0

  # δ_i ≤ a_{i+1}; δ_i ≥ a_{i+1} - b_{i+1} (the map H^{i+1}(A) → H^{i+1}(B)
  # has rank ≥ a_{i+1} - δ_i and at most b_{i+1})
  for i in 0:(d - 1)
    ub[i + 2] = min(ub[i + 2], a_vals[i + 2])
    lb[i + 2] = max(lb[i + 2], a_vals[i + 2] - b_vals[i + 2])
  end

  for _ in 1:(d + 2)
    changed = false

    # c_i ≥ 0 gives δ_{i-1} + δ_i ≥ a_i - b_i; propagate in both directions.
    for i in 0:d
      needed = a_vals[i + 1] - b_vals[i + 1]
      for (j, k) in ((i + 2, i + 1), (i + 1, i + 2))
        new_lb = needed - ub[k]
        if new_lb > lb[j]
          lb[j] = new_lb
          changed = true
        end
      end
    end

    for j in 1:(d + 2)
      lb[j] > ub[j] && (lb[j] = ub[j])
    end

    !changed && break
  end

  (lb, ub)
end

"""
    solve_ses_cohomology(a::Cohomology{BigInt}, b::Cohomology{BigInt})
      -> (Cohomology{BigInt}, Bool)

Given ``H^*(A)`` and ``H^*(B)`` from a short exact sequence
``0 \\to A \\to B \\to C \\to 0``, determine ``H^*(C)`` via the long
exact sequence.

Returns `(H*(C), determined)` where `determined` is `true` when
the long exact sequence uniquely determines all cohomology groups.
When `determined` is `false` the entries are lower bounds and must not
be used as values; use the symbolic solvers instead.
"""
function solve_ses_cohomology(a::Cohomology{BigInt}, b::Cohomology{BigInt})
  d = a.dim_variety
  b.dim_variety == d || throw(
    ArgumentError(
      "solve_ses_cohomology: a and b must have equal dim_variety, got $d and $(b.dim_variety)"
    ),
  )

  a_vals = BigInt[a[i] for i in 0:d]
  b_vals = BigInt[b[i] for i in 0:d]
  (lb, ub) = _ses_delta_bounds(a_vals, b_vals, d)

  # c_i = b_i - a_i + δ_{i-1} + δ_i, using the lower bounds for the δ's.
  entries = BigInt[b_vals[i + 1] - a_vals[i + 1] + lb[i + 1] + lb[i + 2] for i in 0:d]
  (Cohomology{BigInt}(entries, d), lb == ub)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul filtration solver (numeric)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_koszul_filtration(
      koszul_cohos::Vector{Cohomology{BigInt}},
      dim_zero_locus::Int
    ) -> (Cohomology{BigInt}, Bool)

Given ``H^*(X, K_i)`` for the Koszul terms
``K_i = F \\otimes \\wedge^i E^*`` (``i = 0, 1, \\ldots, r``),
compute ``H^*(Z, F|_Z)`` by iterating through the short exact sequences
arising from the Koszul filtration:

```
  C_r = K_r,  0 → C_{j+1} → K_j → C_j → 0   (j = r-1, …, 0)
```

where ``C_0 \\cong F|_Z`` (shifted to the zero locus dimension).

Returns `(H*(F|_Z), determined)`.  When the chain of short exact sequences
leaves some connecting maps open, a symbolic re-solve with the vanishing
``H^k(Z, F|_Z) = 0`` for ``k > \\dim Z`` and the exact Euler characteristic
is attempted before giving up.
"""
function solve_koszul_filtration(
  koszul_cohos::Vector{Cohomology{BigInt}},
  dim_zero_locus::Int,
)
  r = length(koszul_cohos) - 1  # rank of the defining bundle

  current = koszul_cohos[r + 1]  # C_r = K_r
  all_determined = true
  for j in (r - 1):-1:0
    (current, det) = solve_ses_cohomology(current, koszul_cohos[j + 1])
    all_determined = all_determined && det
  end

  result = Cohomology{BigInt}(BigInt[current[i] for i in 0:dim_zero_locus], dim_zero_locus)
  all_determined && return (result, true)

  # Some intermediate SES was undetermined, but the final output may still be
  # pinned down: re-solve symbolically over the full ambient range, then apply
  # the vanishing H^k(Z, F|_Z) = 0 for k > dim Z and the exact Euler
  # characteristic χ(F|_Z) = Σ_i (-1)^i χ(K_i).
  var_counter = Ref(0)
  dim_ambient = koszul_cohos[1].dim_variety
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
#  Short exact sequence solvers (symbolic)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_ses_cohomology_symbolic(a, b, var_counter) -> Cohomology{AffineExpr}

Symbolic version of [`solve_ses_cohomology`](@ref): solve
``0 \\to A \\to B \\to C \\to 0`` for ``H^*(C)``, introducing a fresh
symbolic variable for every connecting-map rank that is not forced.

For numeric inputs (`Cohomology{BigInt}`) the connecting-map bounds are
propagated first and variables are only introduced where a gap remains.
For symbolic inputs each rank ``δ_i`` is zero when ``H^{i+1}(A) = 0``,
equal to ``a_{i+1}`` when ``H^{i+1}(B) = 0`` (the connecting map is then
surjective), and a fresh variable otherwise.

The `var_counter` is advanced for each new variable; pass a shared counter
across multiple calls to keep variable ids globally distinct.
"""
function solve_ses_cohomology_symbolic(
  a::Cohomology{BigInt}, b::Cohomology{BigInt},
  var_counter::Ref{Int},
)
  d = a.dim_variety
  b.dim_variety == d || throw(
    ArgumentError(
      "solve_ses_cohomology_symbolic: a and b must have equal dim_variety, got $d and $(b.dim_variety)"
    ),
  )

  a_vals = BigInt[a[i] for i in 0:d]
  b_vals = BigInt[b[i] for i in 0:d]
  (lb, ub) = _ses_delta_bounds(a_vals, b_vals, d)

  δ = Vector{AffineExpr}(undef, d + 2)
  for j in 1:(d + 2)
    δ[j] = if lb[j] == ub[j]
      AffineExpr(lb[j])
    else
      AffineExpr(lb[j]) + _fresh_variable(var_counter)
    end
  end

  entries = AffineExpr[
    AffineExpr(b_vals[i + 1] - a_vals[i + 1]) + δ[i + 1] + δ[i + 2] for i in 0:d
  ]
  Cohomology{AffineExpr}(entries, d)
end

function solve_ses_cohomology_symbolic(
  a::Cohomology{AffineExpr}, b::Cohomology{AffineExpr},
  var_counter::Ref{Int},
)
  d = a.dim_variety
  b.dim_variety == d || throw(
    ArgumentError(
      "solve_ses_cohomology_symbolic: a and b must have equal dim_variety, got $d and $(b.dim_variety)"
    ),
  )

  # Fully determined inputs: use the numeric bound propagation, it is sharper.
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
      AffineExpr(0)
    elseif is_zero_expr(b[i + 1])
      a[i + 1]  # surjective connecting map
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
  d = b.dim_variety
  b_affine = Cohomology{AffineExpr}(AffineExpr[AffineExpr(b[i]) for i in 0:d], d)
  solve_ses_cohomology_symbolic(a, b_affine, var_counter)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul filtration solver (symbolic)
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

  # C_r = K_r, then 0 → C_{j+1} → K_j → C_j → 0 for j = r-1, …, 0.
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
#  `shortExactSequenceCoker`.
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
Interleave the cohomologies of a short exact sequence
``0 \\to A \\to B \\to C \\to 0`` into its long exact sequence
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

function _numeric_les_cokernel(a::Vector{BigInt}, b::Vector{BigInt})
  d = length(a) - 1
  c, determined = solve_ses_cohomology(Cohomology{BigInt}(a, d), Cohomology{BigInt}(b, d))
  determined || return nothing
  AffineExpr[AffineExpr(c[i]) for i in 0:d]
end

"""
    les_cokernel(a, b, var_counter) -> Vector{AffineExpr}

Given ``H^*(A)`` and ``H^*(B)`` from a short exact sequence
``0 \\to A \\to B \\to C \\to 0``, compute ``H^*(C)``.

Creates a fresh symbolic variable for each ``H^i(C)``, then eliminates as
many as possible using the alternating-sum equations of the long exact
sequence (see `_les_equations`).  Fully determined numeric input is
delegated to the bound-propagation solver first.
"""
function les_cokernel(
  a::Vector{AffineExpr}, b::Vector{AffineExpr},
  var_counter::Ref{Int},
)
  n = length(a)
  length(b) == n || throw(
    ArgumentError("les_cokernel: a and b must have equal length, got $n and $(length(b))")
  )

  all(is_zero_expr, a) && return copy(b)

  if all(is_determined, a) && all(is_determined, b)
    numeric = _numeric_les_cokernel(_determined_bigints(a), _determined_bigints(b))
    numeric !== nothing && return numeric
  end

  c = AffineExpr[_fresh_variable(var_counter) for _ in 1:n]
  for eq in _les_equations(_les_interleave(a, b, c))
    _apply_equation_in_vars!(c, eq)
  end
  c
end

function les_cokernel(a::Vector{AffineExpr}, b::Vector{BigInt}, var_counter::Ref{Int})
  les_cokernel(a, _as_affine(b), var_counter)
end

"""
    les_kernel(b, c, var_counter) -> Vector{AffineExpr}

Given ``H^*(B)`` and ``H^*(C)`` from a short exact sequence
``0 \\to A \\to B \\to C \\to 0``, compute ``H^*(A)``.

Dual to [`les_cokernel`](@ref): the unknowns sit in the ``A``-slots of the
long exact sequence.
"""
function les_kernel(
  b::Vector{AffineExpr}, c::Vector{AffineExpr},
  var_counter::Ref{Int},
)
  n = length(b)
  length(c) == n || throw(
    ArgumentError("les_kernel: b and c must have equal length, got $n and $(length(c))")
  )

  all(is_zero_expr, c) && return copy(b)

  a = AffineExpr[_fresh_variable(var_counter) for _ in 1:n]
  for eq in _les_equations(_les_interleave(a, b, c))
    _apply_equation_in_vars!(a, eq)
  end
  a
end

"""
    long_exact_sequence_cokernel(terms, var_counter) -> Vector{AffineExpr}

Given cohomology of terms `[K_r, K_{r-1}, …, K_0]` (reversed Koszul order),
iteratively apply [`les_cokernel`](@ref) to compute the final cokernel:

```
  C_r = K_r,  0 → C_{j+1} → K_j → C_j → 0  (j = r-1, …, 0)
```
"""
function long_exact_sequence_cokernel(
  terms::Vector{Vector{T}}, var_counter::Ref{Int}
) where {T<:Union{BigInt,AffineExpr}}
  current = _as_affine(terms[1])
  for j in 2:length(terms)
    current = les_cokernel(current, terms[j], var_counter)
  end
  current
end
