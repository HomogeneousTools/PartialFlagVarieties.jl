# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul.jl — Long exact sequence algebra for Koszul resolutions
#
#  Pure algebra module: solves long exact sequences arising from
#  Koszul-type filtrations.  Knows nothing about flag varieties or bundles.
#
#  Given the dimension-valued cohomology of the terms in a short exact
#  sequence or a Koszul filtration, computes the cohomology of the
#  remaining term using the connecting homomorphisms.
# ═══════════════════════════════════════════════════════════════════════════════

export solve_ses_cohomology, solve_koszul_filtration
export AffineExpr, is_determined, is_zero_expr, symbolic_variable
export solve_ses_cohomology_symbolic, solve_koszul_filtration_symbolic

# ═══════════════════════════════════════════════════════════════════════════════
#  Short exact sequence solver
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_ses_cohomology(a::Cohomology{BigInt}, b::Cohomology{BigInt})
      -> (Cohomology{BigInt}, Bool)

Given ``H^*(A)`` and ``H^*(B)`` from a short exact sequence
``0 \\to A \\to B \\to C \\to 0``, determine ``H^*(C)`` via the long
exact sequence.

Returns `(H*(C), determined)` where `determined` is `true` when
the long exact sequence uniquely determines all cohomology groups.

# Algorithm

The long exact sequence reads:
```
⋯ → Hⁱ(A) → Hⁱ(B) → Hⁱ(C) → H^{i+1}(A) → ⋯
```

Denoting ``δ_i = \\mathrm{rank}(H^i(C) \\to H^{i+1}(A))`` (connecting map),
exactness gives:
```
c_i = b_i - a_i + δ_{i-1} + δ_i
```
subject to ``0 \\le δ_i \\le \\min(c_i, a_{i+1})``.

We propagate bounds forward and backward until convergence.
"""
function solve_ses_cohomology(a::Cohomology{BigInt}, b::Cohomology{BigInt})
  d = a.dim_variety
  @assert b.dim_variety == d

  # Initialize bounds on δ_i for i = -1, 0, ..., d
  # δ_{-1} = 0 (no H^{-1}), δ_d = 0 (no H^{d+1}(A))
  lb = zeros(BigInt, d + 2)  # lb[i+2] = lower bound on δ_i, i = -1..d
  ub = fill(BigInt(10)^18, d + 2)

  lb[1] = BigInt(0)  # δ_{-1} = 0
  ub[1] = BigInt(0)
  lb[d + 2] = BigInt(0)  # δ_d = 0
  ub[d + 2] = BigInt(0)

  # δ_i ≤ a_{i+1} for i = 0, ..., d-1
  # δ_i ≥ max(a_{i+1} - b_{i+1}, 0) for i = 0, ..., d-1
  # (rank of connecting map δ_i equals dim ker(H^{i+1}(A)→H^{i+1}(B)),
  #  and that map has rank ≤ min(a_{i+1}, b_{i+1}), so ker ≥ a_{i+1} - b_{i+1})
  for i in 0:(d - 1)
    ub[i + 2] = min(ub[i + 2], a[i + 1])
    new_lb = a[i + 1] - b[i + 1]
    if new_lb > lb[i + 2]
      lb[i + 2] = new_lb
    end
  end

  # Iterate forward-backward passes until convergence
  for _ in 1:(d + 2)
    changed = false

    # Forward pass
    for i in 0:d
      ai = a[i]
      bi = b[i]

      # c_i = b_i - a_i + δ_{i-1} + δ_i ≥ 0
      # δ_i ≥ a_i - b_i - δ_{i-1}
      new_lb = ai - bi - ub[i + 1]
      if new_lb > lb[i + 2]
        lb[i + 2] = new_lb
        changed = true
      end

      # δ_i ≤ b_i - a_i + ... but more useful:
      # c_i ≥ 0 ⟹ δ_{i-1} + δ_i ≥ a_i - b_i
      # c_i ≤ ... but we mainly use δ_i ≤ min(c_i, a_{i+1})
      # c_i = b_i - a_i + δ_{i-1} + δ_i, to have δ_i ≤ c_i we need:
      # δ_i ≤ b_i - a_i + δ_{i-1} + δ_i → always true
      # More useful: from c_i ≤ some bound (none generic)
    end

    # Backward pass
    for i in d:-1:0
      ai = a[i]
      bi = b[i]

      # δ_{i-1} ≥ a_i - b_i - δ_i
      new_lb = ai - bi - ub[i + 2]
      if new_lb > lb[i + 1]
        lb[i + 1] = new_lb
        changed = true
      end
    end

    # Tighten upper bounds from c_i ≥ 0
    for i in 0:d
      ai = a[i]
      bi = b[i]

      # δ_{i-1} + δ_i ≥ a_i - b_i
      total_needed = ai - bi
      if total_needed > 0
        # If ub on one forces lb on the other
        new_lb_prev = total_needed - ub[i + 2]
        if new_lb_prev > lb[i + 1]
          lb[i + 1] = new_lb_prev
          changed = true
        end
        new_lb_curr = total_needed - ub[i + 1]
        if new_lb_curr > lb[i + 2]
          lb[i + 2] = new_lb_curr
          changed = true
        end
      end
    end

    # Clamp lb ≤ ub
    for j in 1:(d + 2)
      if lb[j] > ub[j]
        lb[j] = ub[j]
      end
    end

    !changed && break
  end

  # Check determinacy
  determined = all(lb[j] == ub[j] for j in 1:(d + 2))

  # Compute c_i using the δ values (use lower bounds)
  entries = BigInt[]
  for i in 0:d
    ai = a[i]
    bi = b[i]
    δ_prev = lb[i + 1]
    δ_curr = lb[i + 2]
    ci = bi - ai + δ_prev + δ_curr
    push!(entries, ci)
  end

  (Cohomology{BigInt}(entries, d), determined)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Koszul filtration solver
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_koszul_filtration(
      koszul_cohos::Vector{Cohomology{BigInt}},
      dim_zero_locus::Int
    ) -> (Cohomology{BigInt}, Bool)

Given ``H^*(X, K_i)`` for the Koszul terms
``K_i = F \\otimes \\wedge^i E^*`` (i = 0, 1, \\ldots, r),
compute ``H^*(Z, F|_Z)`` by iterating through the short exact sequences
arising from the Koszul filtration.

The filtration is:
```
  C_r = K_r,  0 → C_{j+1} → K_j → C_j → 0   (j = r-1, …, 0)
```
where ``C_0 \\cong F|_Z`` (shifted to the zero locus dimension).

Returns `(H*(F|_Z), determined)`.
"""
function solve_koszul_filtration(
  koszul_cohos::Vector{Cohomology{BigInt}},
  dim_zero_locus::Int,
)
  r = length(koszul_cohos) - 1  # rank of the defining bundle

  # Start with C_r = K_r
  current = koszul_cohos[r + 1]
  all_determined = true

  # Iterate: 0 → C_{j+1} → K_j → C_j → 0
  for j in (r - 1):-1:0
    (current, det) = solve_ses_cohomology(current, koszul_cohos[j + 1])
    all_determined = all_determined && det
  end

  # current is now C_0 on the ambient X
  # Restrict to the zero locus dimension
  entries = BigInt[]
  for i in 0:dim_zero_locus
    push!(entries, current[i])
  end
  result = Cohomology{BigInt}(entries, dim_zero_locus)

  if all_determined
    return (result, true)
  end

  # Some intermediate SES had undetermined connecting maps, but the final output
  # may still be uniquely determined (uncertain deltas may cancel in later steps).
  # Use the symbolic solver over the FULL ambient dimension and apply vanishing
  # constraints H^k = 0 for k > dim_zero_locus.
  var_counter = Ref(0)
  dim_ambient = koszul_cohos[1].dim_variety
  sym_result_full = solve_koszul_filtration_symbolic(koszul_cohos, dim_ambient, var_counter)

  # Build a matrix over the full ambient range.
  n_full = dim_ambient + 1
  mat = Matrix{AffineExpr}(undef, 1, n_full)
  for k in 0:dim_ambient
    mat[1, k + 1] = sym_result_full[k]
  end

  # Apply vanishing constraints: H^k(Z, F) = 0 for k > dim_zero_locus.
  for k in (dim_zero_locus + 1):dim_ambient
    expr = mat[1, k + 1]
    if !is_determined(expr) || expr.constant != 0
      _apply_equation!(mat, expr - AffineExpr(BigInt(0)))
    end
  end

  # Extract symbolic entries at 0..dim_zero_locus.
  sym_result = Cohomology{AffineExpr}(
    AffineExpr[mat[1, k + 1] for k in 0:dim_zero_locus], dim_zero_locus
  )

  # Apply the Euler characteristic constraint: χ(F|_Z) = Σ_i (-1)^i χ(K_i).
  # χ(K_i) is computed exactly from the ambient Koszul terms.
  chi_exact = _alternating_euler_characteristic(koszul_cohos)

  # Re-use mat from above (restricted to 0..dim_zero_locus)
  mat2 = Matrix{AffineExpr}(undef, 1, dim_zero_locus + 1)
  for k in 0:dim_zero_locus
    mat2[1, k + 1] = sym_result[k]
  end

  # Apply: Σ_k (-1)^k H^k(F|_Z) = chi_exact
  alt_sum = _alternating_sum(mat2, 1, dim_zero_locus)
  _apply_equation!(mat2, alt_sum - AffineExpr(chi_exact))

  # Apply non-negativity propagation: if a symbolic entry equals a determined value,
  # substitute it. Repeat until stable.
  changed = true
  while changed
    changed = false
    for k in 0:dim_zero_locus
      if is_determined(mat2[1, k + 1])
        continue
      end
      # Re-apply χ constraint in case new substitutions freed things
      alt_sum2 = _alternating_sum(mat2, 1, dim_zero_locus)
      changed = _apply_equation!(mat2, alt_sum2 - AffineExpr(chi_exact)) || changed
    end
  end

  sym_determined = all(is_determined(mat2[1, k + 1]) for k in 0:dim_zero_locus)

  if sym_determined
    final_entries = BigInt[mat2[1, k + 1].constant for k in 0:dim_zero_locus]
    final_result = Cohomology{BigInt}(final_entries, dim_zero_locus)
    return (final_result, true)
  end

  (result, false)
end

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

"""Check whether the expression is fully determined (no symbolic variables)."""
is_determined(e::AffineExpr) = isempty(e.coeffs)

"""Check whether the expression is identically zero."""
is_zero_expr(e::AffineExpr) = e.constant == 0 && isempty(e.coeffs)

function _determined_bigints(entries::Vector{AffineExpr})
  BigInt[e.constant for e in entries]
end

function _numeric_les_cokernel(a::Vector{BigInt}, b::Vector{BigInt})
  d = length(a) - 1
  c, determined = solve_ses_cohomology(Cohomology{BigInt}(a, d), Cohomology{BigInt}(b, d))
  determined || return nothing
  AffineExpr[AffineExpr(c[i]) for i in 0:d]
end

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
  c == 1 && return isempty(a.coeffs) ? AffineExpr(a.constant) : AffineExpr(a.constant, copy(a.coeffs))
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
  total = AffineExpr(0)
  add_term = true
  for k in 0:dim
    entry = mat[row, k + 1]
    total = add_term ? total + entry : total - entry
    add_term = !add_term
  end
  total
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
  first_term && print(io, "0")
end

# ─── Cohomology{AffineExpr} display ─────────────────────────────────────────

function Base.show(io::IO, H::Cohomology{AffineExpr})
  parts = String[]
  for i in 0:H.dim_variety
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
#  Symbolic constraint propagation helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""
Substitute ``x_{\\mathrm{var\\_id}} = \\mathrm{replacement}`` in every entry
of the ``\\mathrm{AffineExpr}`` matrix ``M``, in-place.
"""
function _substitute_var!(M::Matrix{AffineExpr}, var_id::Int, replacement::AffineExpr)
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
  end
end

"""
Given a linear equation ``\\mathrm{expr} = 0`` in the symbolic variables,
solve for the variable with the smallest index and substitute the solution
throughout the matrix ``M``.  Returns `true` if a variable was eliminated.

Only eliminates when integer divisibility holds.
"""
function _apply_equation!(M::Matrix{AffineExpr}, expr::AffineExpr)
  isempty(expr.coeffs) && return false
  var_id = minimum(keys(expr.coeffs))
  coeff = expr.coeffs[var_id]

  rest_const = expr.constant
  rest_coeffs = copy(expr.coeffs)
  delete!(rest_coeffs, var_id)

  # Check integer divisibility
  rest_const % coeff == 0 || return false
  all(v % coeff == 0 for (_, v) in rest_coeffs) || return false

  sub_const = -(rest_const ÷ coeff)
  sub_coeffs = Dict{Int,BigInt}(k => -(v ÷ coeff) for (k, v) in rest_coeffs)
  filter!(p -> p.second != 0, sub_coeffs)
  _substitute_var!(M, var_id, AffineExpr(sub_const, sub_coeffs))
  true
end

"""
Apply the constraint ``\\mathrm{hodge}[pi, qi] = \\mathrm{target}`` by
eliminating one symbolic variable.
"""
function _apply_linear_constraint!(
  hodge::Matrix{AffineExpr}, pi::Int, qi::Int, target::BigInt
)
  expr = hodge[pi, qi]
  is_determined(expr) && return false
  _apply_equation!(hodge, expr - AffineExpr(target))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Internal: shared bound propagation
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compute δ bounds for a short exact sequence 0 → A → B → C → 0
given numeric (BigInt) cohomology values for A and B.

Returns `(lb, ub)` — vectors of length `d+2` indexed for δ_{-1}…δ_d.
"""
function _ses_delta_bounds(a_vals::Vector{BigInt}, b_vals::Vector{BigInt}, d::Int)
  lb = zeros(BigInt, d + 2)
  ub = fill(BigInt(10)^18, d + 2)

  lb[1] = BigInt(0);
  ub[1] = BigInt(0)      # δ_{-1} = 0
  lb[d + 2] = BigInt(0);
  ub[d + 2] = BigInt(0)  # δ_d = 0

  for i in 0:(d - 1)
    ub[i + 2] = min(ub[i + 2], a_vals[i + 2])  # δ_i ≤ a_{i+1}
    new_lb = a_vals[i + 2] - b_vals[i + 2]      # δ_i ≥ a_{i+1} - b_{i+1}
    if new_lb > lb[i + 2]
      lb[i + 2] = new_lb
    end
  end

  for _ in 1:(d + 2)
    changed = false

    for i in 0:d
      ai = a_vals[i + 1]
      bi = b_vals[i + 1]
      new_lb = ai - bi - ub[i + 1]
      if new_lb > lb[i + 2]
        lb[i + 2] = new_lb
        changed = true
      end
    end

    for i in d:-1:0
      ai = a_vals[i + 1]
      bi = b_vals[i + 1]
      new_lb = ai - bi - ub[i + 2]
      if new_lb > lb[i + 1]
        lb[i + 1] = new_lb
        changed = true
      end
    end

    for i in 0:d
      ai = a_vals[i + 1]
      bi = b_vals[i + 1]
      total_needed = ai - bi
      if total_needed > 0
        new_lb_prev = total_needed - ub[i + 2]
        if new_lb_prev > lb[i + 1]
          lb[i + 1] = new_lb_prev
          changed = true
        end
        new_lb_curr = total_needed - ub[i + 1]
        if new_lb_curr > lb[i + 2]
          lb[i + 2] = new_lb_curr
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

# ═══════════════════════════════════════════════════════════════════════════════
#  Symbolic short exact sequence solver
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_ses_cohomology_symbolic(
      a::Cohomology{BigInt}, b::Cohomology{BigInt},
      var_counter::Ref{Int}
    ) -> Cohomology{AffineExpr}

Symbolic version of `solve_ses_cohomology`.  When the connecting-map
ranks are not uniquely determined, fresh symbolic variables are
introduced instead of using lower bounds.

The `var_counter` is incremented for each new variable; pass a shared
counter across multiple calls to get globally consistent variable names.
"""
function solve_ses_cohomology_symbolic(
  a::Cohomology{BigInt}, b::Cohomology{BigInt},
  var_counter::Ref{Int},
)
  d = a.dim_variety
  @assert b.dim_variety == d

  a_vals = BigInt[a[i] for i in 0:d]
  b_vals = BigInt[b[i] for i in 0:d]
  (lb, ub) = _ses_delta_bounds(a_vals, b_vals, d)

  # Build symbolic δ values
  δ = Vector{AffineExpr}(undef, d + 2)
  for j in 1:(d + 2)
    if lb[j] == ub[j]
      δ[j] = AffineExpr(lb[j])
    else
      var_counter[] += 1
      δ[j] = AffineExpr(lb[j]) + symbolic_variable(var_counter[])
    end
  end

  # c_i = b_i - a_i + δ_{i-1} + δ_i
  entries = AffineExpr[]
  for i in 0:d
    ci = AffineExpr(b_vals[i + 1] - a_vals[i + 1]) + δ[i + 1] + δ[i + 2]
    push!(entries, ci)
  end

  Cohomology{AffineExpr}(entries, d)
end

"""
    solve_ses_cohomology_symbolic(
      a::Cohomology{AffineExpr}, b::Cohomology{BigInt},
      var_counter::Ref{Int}
    ) -> Cohomology{AffineExpr}

Chained version: A has symbolic entries from a previous step, B is exact.

When A has symbolic entries, cross-bound propagation using constant parts
can be too restrictive (forcing δ_i = 0 when a_{i+1} has non-zero symbolic
parts). Instead, we determine each δ_i directly:
- δ_{-1} = δ_d = 0 (boundary)
- δ_i = 0 if a_{i+1} ≡ 0 (identically zero)
- δ_i = 0 if a_{i+1} is fully determined and equals 0
- δ_i = fresh symbolic variable otherwise

For fully determined inputs (all a_i have no symbolic parts), delegates
to the exact method.
"""
function solve_ses_cohomology_symbolic(
  a::Cohomology{AffineExpr}, b::Cohomology{BigInt},
  var_counter::Ref{Int},
)
  d = a.dim_variety
  @assert b.dim_variety == d

  b_vals = BigInt[b[i] for i in 0:d]

  # If all a values are determined, use exact bound propagation
  if all(is_determined(a[i]) for i in 0:d)
    a_exact = Cohomology{BigInt}(BigInt[a[i].constant for i in 0:d], d)
    return solve_ses_cohomology_symbolic(a_exact, b, var_counter)
  end

  # ── Direct δ assignment with constraint propagation ───────────────────
  δ = Vector{AffineExpr}(undef, d + 2)  # δ[j] for j=1..d+2 corresponds to δ_{j-2}

  # δ_{-1} = 0 and δ_d = 0
  δ[1] = AffineExpr(0)
  δ[d + 2] = AffineExpr(0)

  for i in 0:(d - 1)
    a_next = a[i + 1]
    b_next = b_vals[i + 2]  # b[i+1] in the 0-indexed sense

    if is_zero_expr(a_next)
      # a_{i+1} ≡ 0 ⟹ δ_i = 0
      δ[i + 2] = AffineExpr(0)
    elseif is_determined(a_next) && a_next.constant == 0
      # a_{i+1} is a known zero
      δ[i + 2] = AffineExpr(0)
    elseif b_next == 0
      # b_{i+1} = 0 ⟹ δ_i ≥ a_{i+1} - 0 = a_{i+1} and δ_i ≤ a_{i+1} ⟹ δ_i = a_{i+1}
      δ[i + 2] = a_next
    else
      # General case: a_{i+1} is positive or symbolic with potentially positive values
      var_counter[] += 1
      δ[i + 2] = symbolic_variable(var_counter[])
    end
  end

  # c_i = b_i - a_i + δ_{i-1} + δ_i  (a_i is AffineExpr)
  entries = AffineExpr[]
  for i in 0:d
    ci = AffineExpr(b_vals[i + 1]) - a[i] + δ[i + 1] + δ[i + 2]
    push!(entries, ci)
  end

  Cohomology{AffineExpr}(entries, d)
end

"""
    solve_ses_cohomology_symbolic(
      a::Cohomology{AffineExpr}, b::Cohomology{AffineExpr},
      var_counter::Ref{Int}
    ) -> Cohomology{AffineExpr}

Both A and B have symbolic entries.  Uses direct δ assignment:
δ_i = 0 when a_{i+1} ≡ 0, otherwise a fresh variable.
"""
function solve_ses_cohomology_symbolic(
  a::Cohomology{AffineExpr}, b::Cohomology{AffineExpr},
  var_counter::Ref{Int},
)
  d = a.dim_variety
  @assert b.dim_variety == d

  # If both are fully determined, delegate to exact method
  if all(is_determined(a[i]) for i in 0:d) && all(is_determined(b[i]) for i in 0:d)
    a_exact = Cohomology{BigInt}(BigInt[a[i].constant for i in 0:d], d)
    b_exact = Cohomology{BigInt}(BigInt[b[i].constant for i in 0:d], d)
    return solve_ses_cohomology_symbolic(a_exact, b_exact, var_counter)
  end

  # Direct δ assignment
  δ = Vector{AffineExpr}(undef, d + 2)
  δ[1] = AffineExpr(0)
  δ[d + 2] = AffineExpr(0)

  for i in 0:(d - 1)
    a_next = a[i + 1]
    if is_zero_expr(a_next)
      δ[i + 2] = AffineExpr(0)
    elseif is_determined(a_next) && a_next.constant == 0
      δ[i + 2] = AffineExpr(0)
    else
      var_counter[] += 1
      δ[i + 2] = symbolic_variable(var_counter[])
    end
  end

  entries = AffineExpr[]
  for i in 0:d
    ci = b[i] - a[i] + δ[i + 1] + δ[i + 2]
    push!(entries, ci)
  end

  Cohomology{AffineExpr}(entries, d)
end
# ═══════════════════════════════════════════════════════════════════════════════
#  Symbolic Koszul filtration solver
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_koszul_filtration_symbolic(
      koszul_cohos::Vector{Cohomology{BigInt}},
      dim_zero_locus::Int,
      var_counter::Ref{Int},
    ) -> Cohomology{AffineExpr}

Symbolic version of `solve_koszul_filtration`.  Introduces fresh symbolic
variables for each undetermined connecting-map rank rather than using
lower bounds.
"""
function solve_koszul_filtration_symbolic(
  koszul_cohos::Vector{Cohomology{BigInt}},
  dim_zero_locus::Int,
  var_counter::Ref{Int},
)
  r = length(koszul_cohos) - 1

  # First SES: both inputs are exact BigInt
  if r == 0
    # Trivial: C_0 = K_0, restrict to dim_zero_locus
    entries = AffineExpr[AffineExpr(koszul_cohos[1][i]) for i in 0:dim_zero_locus]
    return Cohomology{AffineExpr}(entries, dim_zero_locus)
  end

  # Start with C_r = K_r (exact)
  # First iteration: 0 → K_r → K_{r-1} → C_{r-1} → 0
  current_sym = solve_ses_cohomology_symbolic(
    koszul_cohos[r + 1], koszul_cohos[r], var_counter
  )

  # Subsequent iterations: A is symbolic, B is exact
  for j in (r - 2):-1:0
    current_sym = solve_ses_cohomology_symbolic(
      current_sym, koszul_cohos[j + 1], var_counter
    )
  end

  # Restrict to zero locus dimension
  entries = AffineExpr[current_sym[i] for i in 0:dim_zero_locus]
  Cohomology{AffineExpr}(entries, dim_zero_locus)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Alternative LES cokernel solver
#
#  Alternative to the δ-based symbolic SES solver.  Introduces a fresh
#  symbolic variable for each *output* entry and derives linear equations
#  from the alternating-sum condition on non-zero segments of the LES.
#  This matches the approach in Macaulay2's `shortExactSequenceCoker`.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _apply_equation_in_vars!(M, eq)

Like `_apply_equation!` but only eliminates variables that appear in `M`.
If the equation involves variables not in `M`, they are treated as
parameters (kept as-is in the solution).
"""
function _apply_equation_in_vars!(M::Matrix{AffineExpr}, eq::AffineExpr)
  isempty(eq.coeffs) && return false

  # Collect variable IDs present in M
  mat_vars = Set{Int}()
  for e in M
    for v in keys(e.coeffs)
      push!(mat_vars, v)
    end
  end

  # Among the equation's variables, pick the smallest that appears in M
  candidates = [v for v in keys(eq.coeffs) if v in mat_vars]
  isempty(candidates) && return false

  var_id = minimum(candidates)
  coeff = eq.coeffs[var_id]

  rest_const = eq.constant
  rest_coeffs = copy(eq.coeffs)
  delete!(rest_coeffs, var_id)

  rest_const % coeff == 0 || return false
  all(v % coeff == 0 for (_, v) in rest_coeffs) || return false

  sub_const = -(rest_const ÷ coeff)
  sub_coeffs = Dict{Int,BigInt}(k => -(v ÷ coeff) for (k, v) in rest_coeffs)
  filter!(p -> p.second != 0, sub_coeffs)
  _substitute_var!(M, var_id, AffineExpr(sub_const, sub_coeffs))
  true
end

"""
    les_cokernel(a, b, var_counter) -> Vector{AffineExpr}

Given `H^*(A)` and `H^*(B)` from a short exact sequence
`0 → A → B → C → 0`, compute `H^*(C)` using the entry-based LES solver.

Creates a fresh symbolic variable for each `H^i(C)`, then derives linear
equations from the alternating-sum condition on non-zero segments of the
interleaved LES `(A₀, B₀, C₀, A₁, B₁, C₁, …)`.  Returns the simplified
`C` entries.
"""
function les_cokernel(
  a::Vector{AffineExpr}, b::Vector{AffineExpr},
  var_counter::Ref{Int},
)
  n = length(a)
  @assert length(b) == n

  all(is_zero_expr, a) && return copy(b)

  if all(is_determined, a) && all(is_determined, b)
    numeric = _numeric_les_cokernel(_determined_bigints(a), _determined_bigints(b))
    numeric !== nothing && return numeric
  end

  # Create fresh variables for c
  c = Vector{AffineExpr}(undef, n)
  for i in 1:n
    var_counter[] += 1
    c[i] = symbolic_variable(var_counter[])
  end

  # Build the interleaved LES: a_1, b_1, c_1, a_2, b_2, c_2, ...
  les = Vector{AffineExpr}(undef, 3n)
  for i in 1:n
    les[3(i - 1) + 1] = a[i]
    les[3(i - 1) + 2] = b[i]
    les[3(i - 1) + 3] = c[i]
  end

  # Split at zeros, form alternating-sum = 0 equations
  equations = AffineExpr[]
  subseq = AffineExpr[]
  for entry in les
    if is_zero_expr(entry)
      if !isempty(subseq)
        alt = _alternating_sum(subseq)
        push!(equations, alt)
        subseq = AffineExpr[]
      end
    else
      push!(subseq, entry)
    end
  end
  if !isempty(subseq)
    alt = _alternating_sum(subseq)
    push!(equations, alt)
  end

  # Solve equations, eliminating only c variables
  mat = reshape(copy(c), 1, n)
  for eq in equations
    _apply_equation_in_vars!(mat, eq)
  end

  AffineExpr[mat[1, i] for i in 1:n]
end

function les_cokernel(
  a::Vector{BigInt}, b::Vector{BigInt},
  var_counter::Ref{Int},
)
  all(iszero, a) && return AffineExpr[AffineExpr(x) for x in b]

  numeric = _numeric_les_cokernel(a, b)
  numeric !== nothing && return numeric

  les_cokernel(
    AffineExpr[AffineExpr(x) for x in a],
    AffineExpr[AffineExpr(x) for x in b],
    var_counter,
  )
end

function les_cokernel(
  a::Vector{AffineExpr}, b::Vector{BigInt},
  var_counter::Ref{Int},
)
  all(is_zero_expr, a) && return AffineExpr[AffineExpr(x) for x in b]

  if all(is_determined, a)
    numeric = _numeric_les_cokernel(_determined_bigints(a), b)
    numeric !== nothing && return numeric
  end

  les_cokernel(a, AffineExpr[AffineExpr(x) for x in b], var_counter)
end

"""
    long_exact_sequence_cokernel(terms, var_counter) -> Vector{AffineExpr}

Given cohomology of terms `[K_r, K_{r-1}, …, K_0]` (reversed Koszul order),
iteratively apply `les_cokernel` to compute the final cokernel.

The filtration is:
  C_r = K_r,  0 → C_{j+1} → K_j → C_j → 0  (j = r-1, …, 0)
"""
function long_exact_sequence_cokernel(
  terms::Vector{Vector{BigInt}},
  var_counter::Ref{Int},
)
  r = length(terms) - 1
  r == 0 && return AffineExpr[AffineExpr(x) for x in terms[1]]

  # First step: both inputs are BigInt
  current = les_cokernel(terms[1], terms[2], var_counter)

  # Subsequent steps: A is AffineExpr, B is BigInt
  for j in 3:length(terms)
    current = les_cokernel(current, terms[j], var_counter)
  end

  current
end

function long_exact_sequence_cokernel(
  terms::Vector{Vector{AffineExpr}},
  var_counter::Ref{Int},
)
  r = length(terms) - 1
  r == 0 && return copy(terms[1])

  current = les_cokernel(terms[1], terms[2], var_counter)
  for j in 3:length(terms)
    current = les_cokernel(current, terms[j], var_counter)
  end

  current
end
