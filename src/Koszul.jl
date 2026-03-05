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
  for i in 0:(d - 1)
    ub[i + 2] = min(ub[i + 2], a[i + 1])
  end

  # Iterate forward-backward passes until convergence
  for _ in 1:d + 2
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

  (Cohomology{BigInt}(entries, dim_zero_locus), all_determined)
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
symbolic_variable(var_id::Int) = AffineExpr(BigInt(0), Dict{Int,BigInt}(var_id => BigInt(1)))

"""Check whether the expression is fully determined (no symbolic variables)."""
is_determined(e::AffineExpr) = isempty(e.coeffs)

"""Check whether the expression is identically zero."""
is_zero_expr(e::AffineExpr) = e.constant == 0 && isempty(e.coeffs)

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
  AffineExpr(a.constant + b.constant, _merge_coeffs(+, a.coeffs, b.coeffs))
end

function Base.:-(a::AffineExpr, b::AffineExpr)
  AffineExpr(a.constant - b.constant, _merge_coeffs(-, a.coeffs, b.coeffs))
end

function Base.:-(a::AffineExpr)
  AffineExpr(-a.constant, Dict{Int,BigInt}(k => -v for (k, v) in a.coeffs))
end

function Base.:*(c::Integer, a::AffineExpr)
  c == 0 && return AffineExpr(0)
  AffineExpr(BigInt(c) * a.constant, Dict{Int,BigInt}(k => BigInt(c) * v for (k, v) in a.coeffs))
end
Base.:*(a::AffineExpr, c::Integer) = c * a

function Base.:(==)(a::AffineExpr, b::AffineExpr)
  a.constant == b.constant && a.coeffs == b.coeffs
end

Base.:+(a::AffineExpr, c::Integer) = a + AffineExpr(c)
Base.:+(c::Integer, a::AffineExpr) = AffineExpr(c) + a
Base.:-(a::AffineExpr, c::Integer) = a - AffineExpr(c)
Base.:-(c::Integer, a::AffineExpr) = AffineExpr(c) - a

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, e::AffineExpr)
  if isempty(e.coeffs)
    print(io, e.constant)
    return
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

  lb[1] = BigInt(0); ub[1] = BigInt(0)      # δ_{-1} = 0
  lb[d + 2] = BigInt(0); ub[d + 2] = BigInt(0)  # δ_d = 0

  for i in 0:(d - 1)
    ub[i + 2] = min(ub[i + 2], a_vals[i + 2])  # δ_i ≤ a_{i+1}
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

  # ── Direct δ assignment (no cross-bound propagation) ──────────────────
  δ = Vector{AffineExpr}(undef, d + 2)  # δ[j] for j=1..d+2 corresponds to δ_{j-2}

  # δ_{-1} = 0 and δ_d = 0
  δ[1] = AffineExpr(0)
  δ[d + 2] = AffineExpr(0)

  for i in 0:(d - 1)
    a_next = a[i + 1]
    if is_zero_expr(a_next)
      # a_{i+1} ≡ 0 ⟹ δ_i = 0
      δ[i + 2] = AffineExpr(0)
    elseif is_determined(a_next) && a_next.constant == 0
      # a_{i+1} is a known zero
      δ[i + 2] = AffineExpr(0)
    else
      # a_{i+1} is positive or symbolic with potentially positive values
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
