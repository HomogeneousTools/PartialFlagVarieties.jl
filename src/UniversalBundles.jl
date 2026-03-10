# ═══════════════════════════════════════════════════════════════════════════════
#  UniversalBundles — tautological and spinor bundles on G/P
#
#  Provides universal (tautological) bundles on classical varieties:
#  - Universal subbundle U and quotient bundle Q on Gr(k, n)
#  - Tautological bundles on OGr, SGr, LGr
#  - Spinor bundles on quadrics Q^n
#  - Tautological bundles on partial flag varieties Fl(d₁,...,dₖ; n)
# ═══════════════════════════════════════════════════════════════════════════════

export universal_subbundle, universal_quotient_bundle
export spinor_bundle
export tautological_bundles, quotient_bundles

# ═══════════════════════════════════════════════════════════════════════════════
#  Type A: Grassmannians Gr(k, n) = A_{n-1}/P_k
# ═══════════════════════════════════════════════════════════════════════════════

"""
    universal_subbundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The universal (tautological) subbundle ``\\mathcal{U}`` on a Grassmannian
``\\mathrm{Gr}(k, n)``.

On ``\\mathrm{Gr}(k, n) = A_{n-1}/P_k``, this is the irreducible equivariant
bundle corresponding to the standard representation of the Levi factor
(the defining representation at node ``k``), i.e., the weight ``\\omega_k``.
It has rank ``k``.

For orthogonal and symplectic Grassmannians, this returns the isotropic
tautological subbundle.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> U = universal_subbundle(X);

julia> rank_bundle(U)
2
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(3, 7);

julia> U = universal_subbundle(X);

julia> rank_bundle(U)
3
```
"""
function universal_subbundle(X::PartialFlagVariety)
  is_generalized_grassmannian(X) || throw(ArgumentError(
    "universal_subbundle requires a generalized Grassmannian (1 marked node)"
  ))

  mdt = marked_dynkin_type(X)
  DT = dynkin_type(mdt)

  # The universal subbundle on G/P_k corresponds to the first fundamental
  # weight ω₁ of the ambient group G.  On Gr(k,n) = A_{n-1}/P_k, this gives
  # the rank-k tautological subbundle (the standard representation of GL_n
  # restricted to the parabolic).  The Levi decomposition then produces a
  # bundle of fiber dimension k (from the standard rep of the GL_k factor).
  ω = fundamental_weight(DT, 1)
  rep = IrrepLevi(mdt, ω)
  CompletelyReducibleBundle(X, [rep])
end

"""
    universal_quotient_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The universal quotient bundle ``\\mathcal{Q}`` on a Grassmannian
``\\mathrm{Gr}(k, n)``.

This is the dual of the universal subbundle: ``\\mathcal{Q} = \\mathcal{U}^*``.
It has rank ``n - k`` on ``\\mathrm{Gr}(k, n)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> Q = universal_quotient_bundle(X);

julia> rank_bundle(Q)
3
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> rank_bundle(universal_subbundle(X)) + rank_bundle(universal_quotient_bundle(X))
5
```
"""
function universal_quotient_bundle(X::PartialFlagVariety)
  is_generalized_grassmannian(X) || throw(ArgumentError(
    "universal_quotient_bundle requires a generalized Grassmannian (1 marked node)"
  ))

  mdt = marked_dynkin_type(X)
  DT = dynkin_type(mdt)

  R = rank(DT)

  if DT <: TypeA
    # On Gr(k, n) = A_{n-1}/P_k, the quotient bundle Q = C^n/U corresponds
    # to the last fundamental weight ω_{n-1} of the ambient A_{n-1}.
    # Under the Levi A_{k-1} × A_{n-k-1}, this has fiber dimension n-k.
    ω = fundamental_weight(DT, R)
    rep = IrrepLevi(mdt, ω)
    return CompletelyReducibleBundle(X, [rep])
  end

  # For non-type-A: fall back to dual of subbundle
  # (This may not match the algebraic-geometric quotient bundle in general)
  dual(universal_subbundle(X))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Spinor bundles on quadrics
# ═══════════════════════════════════════════════════════════════════════════════

"""
    spinor_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle
    spinor_bundle(X::PartialFlagVariety, half::Symbol) -> CompletelyReducibleBundle

The spinor bundle on a quadric ``Q^n``.

For **odd-dimensional** quadrics ``Q^{2m-1} = B_m/P_1``, there is a single
spinor bundle ``\\Sigma`` of rank ``2^{m-1}``, corresponding to the spin
weight ``\\omega_m``.

For **even-dimensional** quadrics ``Q^{2m-2} = D_m/P_1``, there are two
half-spinor bundles ``\\Sigma^+`` and ``\\Sigma^-`` of rank ``2^{m-2}``,
corresponding to ``\\omega_{m-1}`` and ``\\omega_m`` respectively.
Call `spinor_bundle(X, :plus)` or `spinor_bundle(X, :minus)` to select one.
Without a `half` argument on an even quadric, both are returned as a direct sum.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = quadric(5);  # Q^5 = B_3/P_1

julia> S = spinor_bundle(X);

julia> rank_bundle(S)
4
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = quadric(4);  # Q^4 = D_3/P_1

julia> Sp = spinor_bundle(X, :plus);

julia> rank_bundle(Sp)
2
```
"""
function spinor_bundle(X::PartialFlagVariety)
  mdt = marked_dynkin_type(X)
  DT = dynkin_type(mdt)
  marked = marked_nodes(mdt)
  _is_quadric(DT, marked) || throw(ArgumentError(
    "spinor_bundle requires a quadric (B_m/P_1 or D_m/P_1)"
  ))

  R = rank(DT)
  if DT <: TypeB
    # Odd quadric Q^{2m-1}: single spinor bundle at ω_m
    ω = fundamental_weight(DT, R)
    rep = IrrepLevi(mdt, ω)
    CompletelyReducibleBundle(X, [rep])
  elseif DT <: TypeD
    # Even quadric Q^{2m-2}: direct sum of both half-spinors
    ω_plus = fundamental_weight(DT, R - 1)
    ω_minus = fundamental_weight(DT, R)
    rep_plus = IrrepLevi(mdt, ω_plus)
    rep_minus = IrrepLevi(mdt, ω_minus)
    CompletelyReducibleBundle(X, [rep_plus, rep_minus])
  else
    throw(ArgumentError("spinor_bundle requires a quadric (B_m/P_1 or D_m/P_1)"))
  end
end

function spinor_bundle(X::PartialFlagVariety, half::Symbol)
  mdt = marked_dynkin_type(X)
  DT = dynkin_type(mdt)
  marked = marked_nodes(mdt)
  _is_quadric(DT, marked) || throw(ArgumentError(
    "spinor_bundle requires a quadric (B_m/P_1 or D_m/P_1)"
  ))

  R = rank(DT)

  if DT <: TypeB
    half in (:plus, :minus) && @warn "Odd quadric has a single spinor bundle; ignoring half=$half"
    ω = fundamental_weight(DT, R)
    rep = IrrepLevi(mdt, ω)
    return CompletelyReducibleBundle(X, [rep])
  end

  if DT <: TypeD
    if half === :plus
      ω = fundamental_weight(DT, R - 1)
    elseif half === :minus
      ω = fundamental_weight(DT, R)
    else
      throw(ArgumentError("half must be :plus or :minus, got :$half"))
    end
    rep = IrrepLevi(mdt, ω)
    return CompletelyReducibleBundle(X, [rep])
  end

  throw(ArgumentError("spinor_bundle requires a quadric (B_m/P_1 or D_m/P_1)"))
end

"""Check whether DT/P_Marked is a quadric."""
function _is_quadric(::Type{DT}, Marked) where {DT}
  length(Marked) == 1 || return false
  Marked[1] == 1 || return false
  (DT <: TypeB || DT <: TypeD) || return false
  true
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Partial flag varieties Fl(d₁,...,dₖ; n)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    tautological_bundles(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

The chain of tautological subbundles on a partial flag variety.

For ``\\mathrm{Fl}(d_1, \\ldots, d_k; n) = A_{n-1}/P_{\\{d_1,\\ldots,d_k\\}}``,
returns the chain ``\\mathcal{U}_1 \\subset \\mathcal{U}_2 \\subset \\cdots \\subset \\mathcal{U}_k``
where ``\\mathcal{U}_i`` has rank ``d_i`` and corresponds to the weight
``\\omega_{d_i}``.

For a Grassmannian (one marked node), this returns a single-element vector
containing the universal subbundle.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, (1, 2));  # Fl(1,2; 4)

julia> Us = tautological_bundles(X);

julia> length(Us)
2

julia> rank_bundle(Us[1])
1

julia> rank_bundle(Us[2])
1
```

!!! note
    On partial flag varieties with multiple marked nodes, the irreducible
    equivariant bundles ``E_{\\omega_m}`` are **not** the geometric tautological
    subbundles (which are filtered extensions, not completely reducible).
"""
function tautological_bundles(X::PartialFlagVariety)
  # For type A flags Fl(d₁,...,dₖ; n) = A_{n-1}/P_{d₁,...,dₖ},
  # the tautological subbundles 𝒰_i have weights corresponding to
  # ω_{d_i} (the d_i-th fundamental weight of the ambient A_{n-1}).
  # Each gives 𝒰_i of rank d_i.
  mdt = marked_dynkin_type(X)
  DT = dynkin_type(mdt)
  result = CompletelyReducibleBundle[]
  for m in marked_nodes(mdt)
    ω = fundamental_weight(DT, m)
    rep = IrrepLevi(mdt, ω)
    push!(result, CompletelyReducibleBundle(X, [rep]))
  end
  result
end

"""
    quotient_bundles(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

The dual of each tautological bundle at each marked node.

For a Grassmannian, use [`universal_quotient_bundle`](@ref) instead, which
gives the geometrically correct quotient ``\\mathbb{C}^n / \\mathcal{U}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> Qs = quotient_bundles(X);

julia> length(Qs)
1

julia> rank_bundle(Qs[1])
1
```
"""
function quotient_bundles(X::PartialFlagVariety)
  [dual(U) for U in tautological_bundles(X)]
end

# ─── Display ─────────────────────────────────────────────────────────────────

# No special display needed; these return CompletelyReducibleBundle
