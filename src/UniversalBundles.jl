# ═══════════════════════════════════════════════════════════════════════════════
#  UniversalBundles — universal and spinor bundles on G/P
#
#  Provides universal (tautological) bundles on classical varieties:
#  - Universal subbundle U, quotient bundle Q and residual bundle R on isotropic Grassmannians.
#  - Spinor bundles on quadrics Q^n
#  - Tautological bundles on partial flag varieties Fl(d₁,...,dₖ; n)
# ═══════════════════════════════════════════════════════════════════════════════

export universal_subbundle, universal_quotient_bundle, residual_bundle
export spinor_bundle
export is_orthogonal_grassmannian, is_quadric
export tautological_bundles, universal_subbundles

# ═══════════════════════════════════════════════════════════════════════════════
#  Generalized Grassmannians: universal, quotient and residual bundle
# ═══════════════════════════════════════════════════════════════════════════════

# TODO: `universal_subbundle`: need to split documentation for version without and with parameter `i::Int`
"""
    universal_subbundle(X::PartialFlagVariety) -> CompletelyReducibleBundle
    universal_subbundle(X::PartialFlagVariety, i::Int) -> Bundle

On ``\\mathrm{Gr}(k, n) = \\mathrm{A}_{n-1}/P_k``, this is the irreducible equivariant
bundle corresponding to the standard representation of the Levi factor.

For orthogonal and symplectic Grassmannians, this returns the _isotropic_
universal subbundle ``\\mathcal{U}``.

This function is intended for **generalized Grassmannians** (one marked node).
For multi-step flags, use [`tautological_bundles`](@ref) instead; those are
convenient completely reducible building blocks rather than a geometric flag of
subbundles.

The optional argument `i` (default `1`) selects the ``i``-th element from
[`universal_subbundles`](@ref): `universal_subbundle(X, i)` is equivalent to
`universal_subbundles(X)[i]` for `i > 1`. Currently, `i > 1` is only supported
for type ``\\mathrm{A}`` partial flag varieties.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> rank_bundle(universal_subbundle(X))
2
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = OGr(3, 7);

julia> rank_bundle(universal_subbundle(X))
3
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, [1, 2]);  # Fl(1,2; 4)

julia> rank_bundle(universal_subbundle(X, 1))
1

julia> rank_bundle(universal_subbundle(X, 2))
2
```
"""
function universal_subbundle(X::PartialFlagVariety, i::Int=1)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal subbundle")
  )
  DT = dynkin_type(X)
  i == 1 && return dual(CompletelyReducibleBundle(X, fundamental_weight(DT, 1)))
  # TODO: For isotropic Grassmannians (types B, C, D) with i == 2, return U^⊥ as
  # dual(universal_quotient_bundle(X)), reflecting 0 → U → U^⊥ → Res → 0.
  # When Res is zero (Lagrangian/spinor cases), U^⊥ = U, so dual(Q) recovers i == 1.
  k = length(marked_nodes(X))
  (1 ≤ i ≤ k) ||
    throw(ArgumentError("i must be between 1 and length(marked_nodes(X)) = $k"))
  return universal_subbundles(X)[i]
end

"""
    universal_quotient_bundle(X::PartialFlagVariety) -> Union{CompletelyReducibleBundle, FilteredBundle}

TODO: give better explanation, say that this one is completely reducible
On type ``\\mathrm{A}``, this is the equivariant bundle with weight ``\\omega_{n-1}``.
For isotropic Grassmannians (types ``\\mathrm{B}``, ``\\mathrm{C}``, ``\\mathrm{D}``), ``\\mathcal{Q} \\cong \\mathcal{U}^\\vee``
via the bilinear form.
It has rank ``n - k`` on ``\\mathrm{Gr}(k, n)``.

Only the type-``\\mathrm{A}`` case should be read as the literal geometric quotient
``\\mathbb{C}^n / \\mathcal{U}``. For the isotropic cases the implementation
TODO: what is "natural equivariant replacement"? find better explanation!
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

TODO: I don't like this doctest, can we make it more useful? say, compute some H^0?
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> rank_bundle(universal_subbundle(X)) + rank_bundle(universal_quotient_bundle(X))
5
```
"""
function universal_quotient_bundle(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal quotient bundle")
  )
  is_generalized_grassmannian(X) || throw(
    ArgumentError(
      "universal_quotient_bundle requires a generalized Grassmannian (1 marked node)"
    ),
  )

  DT = dynkin_type(X)
  DT <: TypeA && return CompletelyReducibleBundle(X, fundamental_weight(DT, rank(DT)))

  # 0 → U → U^⟂ → R → 0, so Q ≅ dual(U^⟂) (Frassineti–Manivel, arXiv:2605.28712, §1).
  # In Lagrangian / spinor cases R = 0 and U^⟂ = U.
  rank_bundle(residual_bundle(X)) == 0 && return dual(universal_subbundle(X))
  return dual(FilteredBundle(X, [universal_subbundle(X), residual_bundle(X)]))
end

"""
    residual_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

For a generalized Grassmannian of type ``\\mathrm{B}``, ``\\mathrm{C}``, or ``\\mathrm{D}``,
the residual bundle ``\\mathcal{R} = \\mathcal{U}^\\perp / \\mathcal{U}`` is the
completely reducible bundle fitting into the short exact sequence
``0 \\to \\mathcal{U} \\to \\mathcal{U}^\\perp \\to \\mathcal{R} \\to 0``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = OGr(3, 9);

julia> R = residual_bundle(X);

julia> rank_bundle(R)
3
```

TODO: don't write somewhat silly tests like this
TODO: is there maybe a cohomology computation we can verify? e.g., the Ext^1 of the pieces?
```jldoctest
julia> using PartialFlagVarieties

julia> X = SGr(3, 8);

julia> rank_bundle(residual_bundle(X)) + rank_bundle(universal_subbundle(X)) == rank_bundle(universal_quotient_bundle(X))
true
```

"""
function residual_bundle(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined residual bundle")
  )
  is_generalized_grassmannian(X) || throw(
    ArgumentError(
      "residual_bundle requires a generalized Grassmannian (1 marked node)"
    ),
  )
  DT = dynkin_type(X)
  DT <: TypeA &&
    throw(ArgumentError("type A has no well-defined residual bundle"))

  marked = marked_nodes(X)[1]
  R = rank(DT)

  marked == 1 && throw(
    ArgumentError("not implemented for marked node 1")
  )

  # Weight formulas: Frassineti–Manivel, arXiv:2605.28712, §1.
  if DT <: TypeB
    if marked == R
      ω = WeightLatticeElem(DT)
    elseif marked == R - 1
      ω = 2 * fundamental_weight(DT, R) - fundamental_weight(DT, R - 1)
    else
      ω = fundamental_weight(DT, marked + 1) - fundamental_weight(DT, marked)
    end
  elseif DT <: TypeC
    if marked == R
      return zero_bundle(X)
    else
      ω = fundamental_weight(DT, marked + 1) - fundamental_weight(DT, marked)
    end
  elseif DT <: TypeD
    if marked in (R, R - 1)
      return zero_bundle(X)
    elseif marked == R - 2
      ω =
        fundamental_weight(DT, R) + fundamental_weight(DT, R - 1) -
        fundamental_weight(DT, R - 2)
    else
      ω = fundamental_weight(DT, marked + 1) - fundamental_weight(DT, marked)
    end
  end
  return CompletelyReducibleBundle(X, ω)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Spinor bundles on orthogonal Grassmannians
# ═══════════════════════════════════════════════════════════════════════════════

"""
    spinor_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle
    spinor_bundle(X::PartialFlagVariety, half::Symbol) -> CompletelyReducibleBundle

The spinor bundle on an orthogonal Grassmannian ``\\mathrm{OGr}(k, m)``.

For type ``\\mathrm{B}_n/P_k`` (``\\mathrm{OGr}(k, 2n+1)``) there is a single
spinor bundle ``\\mathcal{S}`` corresponding to ``\\omega_n``.

For type ``\\mathrm{D}_n/P_k`` (``\\mathrm{OGr}(k, 2n)``) there are two
half-spinor bundles ``\\mathcal{S}^+`` and ``\\mathcal{S}^-`` corresponding to
``\\omega_{n-1}`` and ``\\omega_n`` respectively. Call
`spinor_bundle(X, :plus)` or `spinor_bundle(X, :minus)` to select one. Without
a `half` argument on type ``\\mathrm{D}``, both are returned as a direct sum.

For ``k \\le n-2`` (type ``\\mathrm{D}_n``) or ``k \\le n-1`` (type
``\\mathrm{B}_n``), the spinor bundles are the irreducible homogeneous bundles
whose fibre at ``[U]`` is a (half-)spin representation of the residual quadratic
space ``U^\\perp / U``. On the spinor varieties ``\\mathrm{OGr}(n, 2n)``, one
half-spinor bundle is the hyperplane line bundle of the spinor embedding and
the other is a twist of the tautological rank-``n`` bundle.

This constructor additionally accepts the two-marked variety
``\\mathrm{D}_n/P_{n-1,\\,n} = \\mathrm{OGr}(n-1, 2n)`` — the ``(n-1)``-isotropic
Grassmannian, of Picard rank 2. There ``\\mathcal{S}^{\\pm}`` are the
pull-backs of the hyperplane line bundles from ``\\mathrm{OGr}(n,2n)_{\\pm}``
along the embedding
``\\mathrm{OGr}(n-1, 2n) \\hookrightarrow \\mathrm{OGr}(n, 2n)_+ \\times
\\mathrm{OGr}(n, 2n)_-`` that sends an ``(n-1)``-isotropic subspace to the
unique pair of maximal isotropics (one from each family) containing it.

# References

- Ottaviani, "Spinor bundles on quadrics", *Trans. Amer. Math. Soc.* **307**
  (1988), 301–316 — the case ``k = 1`` (quadrics).
- Manivel, "On spinor varieties and their secants", *SIGMA* **5** (2009), 078
  — the spinor varieties ``\\mathrm{OGr}(n, 2n)_\\pm``.
- Frassineti–Manivel, "Spinorial Fano manifolds", arXiv:2605.28712, §1
  ("Spin representations and spin bundles") — general ``\\mathrm{OGr}(k, m)``
  including the ``(n-1)``-isotropic case
  ``\\mathrm{D}_n/P_{n-1, n} = \\mathrm{OGr}(n-1, 2n)`` (Picard rank 2,
  the two spinor bundles pulled back from the maximal spinor varieties).
- Fulton–Harris, *Representation Theory: A First Course*, GTM 129
  (Springer, 1991), §20 — the underlying spin representations of
  ``\\mathrm{Spin}_m``.

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

TODO: use some cohomology doctests too
```jldoctest
julia> using PartialFlagVarieties

julia> X = OGr(3, 10);  # D_5/P_3

julia> rank_bundle(spinor_bundle(X, :plus))
2
```
"""
function spinor_bundle(X::PartialFlagVariety)
  is_orthogonal_grassmannian(X) || throw(
    ArgumentError(
      "spinor_bundle requires an orthogonal Grassmannian (B_n/P_k, D_n/P_k, or D_n/P_{n-1,n})"
    ),
  )
  DT = dynkin_type(X)
  R = rank(DT)
  # OGr(k, 2n+1): single spinor bundle at ω_n.
  DT <: TypeB && return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  # DT <: TypeD by is_orthogonal_grassmannian; OGr(k, 2n) carries both half-spinors.
  return CompletelyReducibleBundle(
    X, [fundamental_weight(DT, R - 1), fundamental_weight(DT, R)]
  )
end

function spinor_bundle(X::PartialFlagVariety, half::Symbol)
  is_orthogonal_grassmannian(X) || throw(
    ArgumentError(
      "spinor_bundle requires an orthogonal Grassmannian (B_n/P_k, D_n/P_k, or D_n/P_{n-1,n})"
    ),
  )
  DT = dynkin_type(X)
  R = rank(DT)

  if DT <: TypeB
    half in (:plus, :minus) &&
      @warn "Type B has a single spinor bundle; ignoring half=$half"
    return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  end

  # DT <: TypeD by is_orthogonal_grassmannian.
  half === :plus && return CompletelyReducibleBundle(X, fundamental_weight(DT, R - 1))
  half === :minus && return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  throw(ArgumentError("half must be :plus or :minus, got :$half"))
end

# TODO: put these checks near `    is_generalized_grassmannian(X::PartialFlagVariety) -> Bool
` in src/PartialFlagVarieties.jl
"""
    is_orthogonal_grassmannian(X::PartialFlagVariety) -> Bool

Return `true` if `X` is an orthogonal Grassmannian, i.e. one of:

- a single-marked variety of type ``\\mathrm{B}_n`` or ``\\mathrm{D}_n``,
  meaning ``\\mathrm{OGr}(k, 2n+1)`` or ``\\mathrm{OGr}(k, 2n)``;
- the two-marked ``\\mathrm{D}_n/P_{n-1, n} = \\mathrm{OGr}(n-1, 2n)``,
  the ``(n-1)``-isotropic Grassmannian of Picard rank 2 (see Frassineti–Manivel,
  arXiv:2605.28712, §1).

This is exactly the class of varieties on which [`spinor_bundle`](@ref) is defined.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_orthogonal_grassmannian(OGr(3, 10))
true

# TODO: is this really an orthogonal Grassmannian?! in the spinor bundle code: maybe add this as a special allowed case, rather than calling it an orthogonal Grassmannian
julia> is_orthogonal_grassmannian(partial_flag_variety(TypeD{4}, (3, 4)))
true

julia> is_orthogonal_grassmannian(Gr(2, 5))
false
```
"""
function is_orthogonal_grassmannian(X::PartialFlagVariety)
  DT = dynkin_type(X)
  (DT <: TypeB || DT <: TypeD) || return false
  marked = marked_nodes(X)
  length(marked) == 1 && return true
  if DT <: TypeD
    R = rank(DT)
    length(marked) == 2 && Set(marked) == Set((R - 1, R)) && return true
  end
  return false
end

"""
    is_quadric(X::PartialFlagVariety) -> Bool

Return `true` if `X` is a smooth quadric hypersurface ``Q^n``, i.e.
``\\mathrm{B}_m/P_1`` (odd ``n = 2m-1``) or ``\\mathrm{D}_m/P_1`` (even ``n = 2m-2``).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_quadric(quadric(4))
true

julia> is_quadric(OGr(2, 10))
false
```
"""
function is_quadric(X::PartialFlagVariety)
  DT = dynkin_type(X)
  (DT <: TypeB || DT <: TypeD) || return false
  marked_nodes(X) == (1,)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Partial flag varieties Fl(d₁,...,dₖ; n)
# ═══════════════════════════════════════════════════════════════════════════════
# TODO: Extend to isotropic types. For a generalized Grassmannian of type B, C, or D
# the building blocks are universal_subbundle(X) and residual_bundle(X). For multi-step
# isotropic flags, compute the graded pieces U_i/U_{i-1} via weight differences
# analogously to the type A case below.
"""
    tautological_bundles(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

Compute the tautological bundles on a partial flag variety ``Fl(d_1, \\ldots, d_k; n)``.

For type ``\\mathrm{A}`` partial flags, returns a vector of completely reducible bundles, where
each element corresponds to a tautological subbundle.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4,[1,2]);  # Flag variety Fl(1,2; 4)

# TODO: weird choice of variable name; why not U?
julia> τ = tautological_bundles(X);

julia> length(τ)
2
```
"""
function tautological_bundles(X::PartialFlagVariety)
  if is_exceptional_type(X)
    throw(
      ArgumentError("exceptional types do not have well-defined tautological bundles")
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
          X, -fundamental_weight(DT, marked[i]) + fundamental_weight(DT, marked[i] - 1)
        ),
      )
    end
    return bundles
  else
    throw(
      ArgumentError(
        "tautological_bundles not implemented for non-type A partial flag varieties"
      ),
    ) #TODO: For isotropic types B, C, D, the graded pieces of the tautological filtration are
    # universal_subbundle(X) and residual_bundle(X) for Grassmannians, and more generally
    # the U_i/U_{i-1} computed from the weight differences at each marked node.
  end
end

# TODO: shouldn't this be closer to universal_subbundle?
"""
    universal_subbundles(X::PartialFlagVariety) -> Vector{Bundle}

Compute the universal subbundles on a partial flag variety as filtered bundles.

For type ``\\mathrm{A}`` partial flags ``Fl(d_1, \\ldots, d_k; n)``, returns a vector of filtered
bundles where each element represents a nested subbundle. The first element is the
first universal bundle (rank ``d_1``), and the ``i``-th element is a filtered bundle
corresponding with the ``i``-th geometric universal subbundle.

These form the complete flag filtration:
``0 \\subset U_1 \\subset U_2 \\subset \\cdots \\subset U_k``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4,[1,2]);  # Flag variety Fl(1,2; 4)

julia> subs = universal_subbundles(X);

julia> rank_bundle(subs[1])
1

julia> rank_bundle(subs[2])
2
```
"""
function universal_subbundles(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have well-defined universal subbundles")
  )
  DT = dynkin_type(X)
  !(DT <: TypeA) && throw(
    ArgumentError(
      "universal_subbundles currently only implemented for type A partial flag varieties"
    ),
  )#TODO: For isotropic Grassmannians (types B, C, D), the filtration steps are
  # [U] and [U^⊥ = FilteredBundle(X, [universal_subbundle(X), residual_bundle(X)])].
  # For isotropic partial flags, extend similarly to the type A implementation above.
  marked = marked_nodes(X)
  bundles = Vector{Bundle}()
  # TODO: don't use this variable name
  τ = tautological_bundles(X)
  push!(bundles, τ[1])
  for i in 2:length(marked)
    push!(bundles, FilteredBundle(X, τ[1:i]))
  end
  return bundles
end
