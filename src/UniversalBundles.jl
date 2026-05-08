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
export tautological_bundles, universal_subbundles

# ═══════════════════════════════════════════════════════════════════════════════
#  Generalized Grassmannians: universal, quotient and residual bundle
# ═══════════════════════════════════════════════════════════════════════════════
"""
    universal_subbundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

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

julia> X = OGr(3, 7);

julia> U = universal_subbundle(X);

julia> rank_bundle(U)
3
```
"""
function universal_subbundle(X::PartialFlagVariety, n::Int=1)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal subbundle")
  )
  DT = dynkin_type(X)
  R = rank(DT)
  if n==1
    return dual(CompletelyReducibleBundle(X, fundamental_weight(DT, 1)))
  else #TODO: If n==2 and X is isotropic Grassmannian, return the dual of the quotient bundle.
    (n < 1 || n > R) &&
      throw(ArgumentError("n must be between 1 and the rank of the variety"))
    return universal_subbundles(X)[n]
  end
end

"""
    universal_quotient_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle || FilteredBundle

On type A, this is the equivariant bundle with weight ``\\omega_{n-1}``.
For isotropic Grassmannians (types B, C, D), ``\\mathcal{Q} \\cong \\mathcal{U}^\\vee``
via the bilinear form.
It has rank ``n - k`` on ``\\mathrm{Gr}(k, n)``.

Only the type-``A`` case should be read as the literal geometric quotient
``\\mathbb{C}^n / \\mathcal{U}``. For the isotropic cases the implementation
returns the dual of the tautological bundle as the natural equivariant
replacement.

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
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal quotient bundle")
  )

  DT = dynkin_type(X)
  R = rank(DT)

  if DT <: TypeA
    return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  else
    # U = universal_subbundle(X), R = residual_bundle(X). In this case the quotient Q^* = U^⟂ bundle fits into the s.e.s
    # 0 -> U -> U^⟂ -> R -> 0
    # So, the code returns the dual of a filtered bundle with subbundle U and quotient R.
    if rank_bundle(residual_bundle(X)) == 0
      return dual(universal_subbundle(X))
    else
      return dual(FilteredBundle(X, [universal_subbundle(X), residual_bundle(X)]))
    end
  end
end
"""
    residual_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

For types B, C, D, the residual bundle is the completely reducible bundle corresponding to U^⟂ / U.

# Examples:
```jldoctest
julia> using PartialFlagVarieties

julia> X = OGr(3, 9);

julia> R = residual_bundle(X);

julia> rank_bundle(R)
3
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = SGr(3, 8);

julia> rank_bundle(residual_bundle(X)) + rank_bundle(universal_subbundle(X)) == rank_bundle(universal_quotient_bundle(X))
true
```

"""
function residual_bundle(X)
  is_generalized_grassmannian(X) || throw(
    ArgumentError(
      "residual_bundle requires a generalized Grassmannian (1 marked node)"
    ),
  )
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined residual bundle")
  )

  DT = dynkin_type(X)
  marked = marked_nodes(X)[1]
  R = rank(DT)

  marked == 1 && throw(
    ArgumentError("not implemented for marked node 1")
  )
  if DT <: TypeA
    throw(ArgumentError("type A do not have a well-defined residual bundle"))
  elseif DT <: TypeB
    if marked == R
      ω = WeightLatticeElem(DT)
    elseif marked == R - 1
      ω = 2 * fundamental_weight(DT, R) - fundamental_weight(DT, R - 1)
    else
      ω = fundamental_weight(DT, marked+1) - fundamental_weight(DT, marked)
    end
  elseif DT <: TypeC
    if marked == R
      return zero_bundle(X)
    else
      ω = fundamental_weight(DT, marked+1) - fundamental_weight(DT, marked)
    end
  elseif DT <: TypeD
    if marked in (R, R-1)
      return zero_bundle(X)
    elseif marked == R - 2
      ω =
        fundamental_weight(DT, R) + fundamental_weight(DT, R - 1) -
        fundamental_weight(DT, R - 2)
    else
      ω = fundamental_weight(DT, marked+1) - fundamental_weight(DT, marked)
    end
  end
  return CompletelyReducibleBundle(X, ω)
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
#TODO Change name of tautological_bundles?
"""
    tautological_bundles(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

Compute the tautological bundles on a partial flag variety ``Fl(d_1, \\ldots, d_k; n)``.

For **type A** partial flags, returns a vector of completely reducible bundles, where
each element corresponds to a tautological subbundle.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Fl(1, 2; 4);  # Flag variety Fl(1,2; 4)

julia> τ = tautological_bundles(X);

julia> length(τ)
2
```
"""
function tautological_bundles(X::PartialFlagVariety)
  if is_exceptional_type(X)
    throw(
      ArgumentError("exceptional types do not have a well-defined tautological bundles")
    )
  end
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  if DT <: TypeA
    bundles = Vector{CompletelyReducibleBundle}()
    push!(bundles, dual(CompletelyReducibleBundle(X, fundamental_weight(DT, 1))))
    for i in 2:length(marked)
      push!(
        bundles,
        CompletelyReducibleBundle(
          X, -fundamental_weight(DT, marked[i]) + fundamental_weight(DT, marked[i]-1)
        ),
      )
    end
    return bundles
  else
    throw(
      ArgumentError(
        "tautological_bundles not implemented for non-type A partial flag varieties"
      ),
    ) #TODO: Implement for other types.
  end
end

"""
    universal_subbundles(X::PartialFlagVariety) -> Vector{Bundle}

Compute the universal subbundles on a partial flag variety as filtered bundles.

For **type A** partial flags ``Fl(d_1, \\ldots, d_k; n)``, returns a vector of filtered
bundles where each element represents a nested subbundle. The first element is the
first universal bundle (rank ``d_1``), and the ``i``-th element is a filtered bundle
corresponding with the ``i``-th geometric universal subundle.

These form the complete flag filtration:
``0 \\subset U_1 \\subset U_2 \\subset \\cdots \\subset U_k``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Fl(1, 2; 4);  # Flag variety Fl(1,2; 4)

julia> subs = universal_subbundles(X);

julia> rank_bundle(subs[1])
1

julia> rank_bundle(subs[2])
2
```
"""
function universal_subbundles(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal subbundles")
  )
  DT = dynkin_type(X)
  !(DT <: TypeA) && throw(
    ArgumentError(
      "universal_subbundles currently only implemented for type A partial flag varieties"
    ),
  )
  marked = marked_nodes(X)
  bundles = Vector{Bundle}()
  τ = tautological_bundles(X)
  push!(bundles, τ[1])
  for i in 2:length(marked)
    push!(bundles, FilteredBundle(X, τ[1:i]))
  end
  return bundles
end
