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
⋯ → H^i(A) → H^i(B) → H^i(C) →  H^{i+1}(A) → ⋯
```

Denoting ``δ_i = \\mathrm{rank}(H^i(C) \\to H^{i+1}(A))`` (connecting map),
exactness gives:
```
c_i = b_i - a_i + δ_{i-1} + δ_i
```
subject to ``0 \\le δ_i \\le \\min(c_i, a_{i+1})``.

We propagate bounds forward and backward to determine all ``δ_i``.
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

  # Forward pass: propagate constraints
  for i in 0:d
    ai = a[i]
    bi = b[i]

    # From c_i ≥ 0: δ_{i-1} + δ_i ≥ a_i - b_i
    # From δ_i ≤ a_{i+1}: upper bound
    if i < d
      ub[i + 2] = min(ub[i + 2], a[i + 1])
    else
      ub[i + 2] = BigInt(0)  # δ_d = 0
    end

    # c_i = b_i - a_i + δ_{i-1} + δ_i ≥ 0
    # δ_i ≥ a_i - b_i - δ_{i-1}
    # With known δ_{i-1} bounds:
    min_delta_prev = lb[i + 1]
    needed = ai - bi - min_delta_prev
    if needed > 0
      lb[i + 2] = max(lb[i + 2], needed)
    end
  end

  # Backward pass: from δ_d = 0
  for i in d:-1:0
    ai = a[i]
    bi = b[i]

    # c_i = b_i - a_i + δ_{i-1} + δ_i ≥ 0
    # δ_{i-1} ≥ a_i - b_i - δ_i
    max_delta_i = ub[i + 2]
    needed = ai - bi - max_delta_i
    if needed > 0
      lb[i + 1] = max(lb[i + 1], needed)
    end
  end

  # Clamp lb ≤ ub (forward pass can overshoot when δ_d = 0 conflicts)
  for i in 1:(d + 2)
    lb[i] = min(lb[i], ub[i])
  end

  # Check determinacy
  determined = all(lb[i] == ub[i] for i in 1:(d + 2))

  if !determined
    # Try to resolve by using c_i ≥ 0 more aggressively
    # Additional forward-backward iterations
    for _ in 1:3
      for i in 0:d
        ai = a[i]
        bi = b[i]
        # c_i = b_i - a_i + δ_{i-1} + δ_i
        # Must have c_i ≥ 0, so δ_{i-1} + δ_i ≥ a_i - b_i
        total_needed = ai - bi
        if total_needed > 0
          # If ub on one forces lb on the other
          if ub[i + 2] < total_needed - lb[i + 1]
            lb[i + 1] = max(lb[i + 1], total_needed - ub[i + 2])
          end
          if ub[i + 1] < total_needed - lb[i + 2]
            lb[i + 2] = max(lb[i + 2], total_needed - ub[i + 1])
          end
        end
      end
      # Clamp lb ≤ ub after each iteration
      for j in 1:(d + 2)
        lb[j] = min(lb[j], ub[j])
      end
      determined = all(lb[j] == ub[j] for j in 1:(d + 2))
      determined && break
    end
  end

  # Compute c_i using the δ values (use lower bounds when determined)
  entries = BigInt[]
  for i in 0:d
    ai = a[i]
    bi = b[i]
    δ_prev = lb[i + 1]
    δ_curr = lb[i + 2]
    ci = bi - ai + δ_prev + δ_curr
    push!(entries, ci)
  end

  return (Cohomology{BigInt}(entries, d), determined)
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

  return (Cohomology{BigInt}(entries, dim_zero_locus), all_determined)
end
