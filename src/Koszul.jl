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
