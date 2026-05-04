# ═══════════════════════════════════════════════════════════════════════════════
#  UniversalBundles — tautological and spinor bundles on G/P
#
#  Provides universal (tautological) bundles on classical varieties:
#  - Universal subbundle U and quotient bundle Q on Gr(k, n)
#  - Tautological bundles on OGr, SGr, LGr
#  - Spinor bundles on quadrics Q^n
#  - Tautological bundles on partial flag varieties Fl(d₁,...,dₖ; n)
# ═══════════════════════════════════════════════════════════════════════════════

export universal_subbundle, universal_quotient_bundle, residual_bundle
export spinor_bundle
export tautological_bundles, quotient_bundles

# ═══════════════════════════════════════════════════════════════════════════════
#  Generalized Grassmannians: universal, quotient and residual bundle
# ═══════════════════════════════════════════════════════════════════════════════

"""
    universal_subbundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The universal (tautological) subbundle ``\\mathcal{U}`` on a Grassmannian
``\\mathrm{Gr}(k, n)``.

On ``\\mathrm{Gr}(k, n) = A_{n-1}/P_k``, this is the irreducible equivariant
bundle corresponding to the standard representation of the Levi factor.

For orthogonal and symplectic Grassmannians, this returns the isotropic
tautological subbundle.

This function is intended for **generalized Grassmannians** (one marked node).
For multi-step flags, use [`tautological_bundles`](@ref) instead; those are
convenient completely reducible building blocks rather than a geometric flag of
subbundles.

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
  is_generalized_grassmannian(X) || throw(
    ArgumentError(
      "universal_subbundle requires a generalized Grassmannian (1 marked node)"
    ),
  )

  mdt = marked_dynkin_type(X)
  is_exceptional(mdt) && throw(
    ArgumentError("exceptional types do not have a well-defined universal subbundle")
  )
  # The equivariant bundle with weight ω₁ is the dual U^∨ of the
  # tautological subbundle (it has global sections = the standard
  # representation, whereas U has none).  The universal subbundle is
  # therefore the dual of this bundle.
  dual(CompletelyReducibleBundle(X, fundamental_weight(dynkin_type(mdt), 1)))
end

"""
    universal_quotient_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The universal quotient bundle ``\\mathcal{Q}`` on a Grassmannian
``\\mathrm{Gr}(k, n)``.

On type A, this is the equivariant bundle with weight ``\\omega_{n-1}``.
For isotropic Grassmannians (types B, C, D), ``\\mathcal{Q} \\cong \\mathcal{U}^\\vee``
via the bilinear form.
It has rank ``n - k`` on ``\\mathrm{Gr}(k, n)``.

Only the type-``A`` case should be read as the literal geometric quotient
``\\mathbb{C}^n / \\mathcal{U}``. For the isotropic cases the implementation
returns the dual of the tautological bundle as the natural equivariant
replacement.

TODO: implement in other types

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
  is_generalized_grassmannian(X) || throw(
    ArgumentError(
      "universal_quotient_bundle requires a generalized Grassmannian (1 marked node)"
    ),
  )

  mdt = marked_dynkin_type(X)
  is_exceptional(mdt) && throw(
    ArgumentError("exceptional types do not have a well-defined universal quotient bundle")
  )

  DT = dynkin_type(mdt)
  R = rank(DT)

  if DT <: TypeA
    # On Gr(k, n) = A_{n-1}/P_k, the quotient bundle Q = C^n/U corresponds
    # to the last fundamental weight ω_{n-1} of the ambient A_{n-1}.
    # Under the Levi A_{k-1} × A_{n-k-1}, this has fiber dimension n-k.
    return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  end

  throw(
    ArgumentError(
      "universal_quotient_bundle is currently only implemented for type A Grassmannians"
    )
  )
end

function residual_bundle(X)
  is_generalized_grassmannian(X) || throw(
    ArgumentError(
      "residual_bundle requires a generalized Grassmannian (1 marked node)"
    ),
  )

  mdt = marked_dynkin_type(X)
  is_exceptional(mdt) && throw(
    ArgumentError("exceptional types do not have a well-defined residual bundle")
  )

  DT = dynkin_type(mdt)
  marked = marked_nodes(mdt)
  R = rank(DT)

  if DT <: TypeB && !(marked[1] in (R, R - 1))
    return CompletelyReducibleBundle(
      X,
      fundamental_weight(DT, marked[1] + 1) - fundamental_weight(DT, marked[1]),
    )
  end

    if DT <: TypeC && marked[1] != R
    return CompletelyReducibleBundle(
      X,
      fundamental_weight(DT, marked[1] + 1) - fundamental_weight(DT, marked[1]),
    )
  end
  throw(ArgumentError("residual_bundle: unsupported type or marked node"))
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

This function is only defined on quadrics, i.e. on ``B_m/P_1`` or ``D_m/P_1``.

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
    half in (:plus, :minus) &&
      @warn "Odd quadric has a single spinor bundle; ignoring half=$half"
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
returns one completely reducible bundle attached to each marked node
``d_i``, corresponding to the ambient fundamental weight ``\\omega_{d_i}``.

For a Grassmannian (one marked node), this returns a single-element vector
containing the universal subbundle and agrees with the geometric tautological
bundle.

For multi-step flags, these are **not** the literal geometric tautological
subbundles ``\\mathcal{U}_i`` inside a nested flag. The geometric bundles are
filtered objects, while this function returns the completely reducible pieces
that are most useful for representation-theoretic computations.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, (1, 2));  # Fl(1,2; 4)

julia> Us = tautological_bundles(X);

julia> length(Us)
2
```

!!! note
    On partial flag varieties with multiple marked nodes, the irreducible
    equivariant bundles ``E_{\\omega_m}`` are **not** the geometric tautological
    subbundles, which are filtered extensions rather than completely reducible
    bundles.
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

For multi-step flags this should again be interpreted as a collection of
convenient completely reducible bundles attached to the marked nodes, not as
the literal successive quotients of a tautological flag.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> Qs = quotient_bundles(X);

julia> length(Qs)
1
```
"""
function quotient_bundles(X::PartialFlagVariety)
  [dual(U) for U in tautological_bundles(X)]
end

# ─── Display ─────────────────────────────────────────────────────────────────

# No special display needed; these return CompletelyReducibleBundle
