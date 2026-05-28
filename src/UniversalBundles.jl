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
export tautological_bundles, universal_subbundles

# ═══════════════════════════════════════════════════════════════════════════════
#  Generalized Grassmannians: universal, quotient and residual bundle
# ═══════════════════════════════════════════════════════════════════════════════

"""
    universal_subbundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

On a generalized Grassmannian ``\\mathrm{Gr}(k, n) = \\mathrm{A}_{n-1}/P_k``,
the irreducible equivariant bundle ``\\mathcal{U}`` corresponding to the
standard representation of the Levi factor — geometrically the tautological
rank-``k`` subbundle of the trivial bundle ``\\mathbb{C}^n``.

For orthogonal and symplectic Grassmannians, returns the _isotropic_ universal
subbundle ``\\mathcal{U}``.

For multi-step partial flag varieties, see [`universal_subbundle(X, i)`](@ref)
and [`universal_subbundles`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(universal_subbundle(Gr(2, 5)))
2
```

```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(universal_subbundle(OGr(3, 7)))
3
```

``\\mathcal{U}`` sits inside the trivial bundle, so ``H^0(\\mathcal{U}) = 0``,
while ``H^0(\\mathcal{U}^\\vee)`` recovers the standard representation of the
structure group:
```jldoctest
julia> using PartialFlagVarieties

julia> U = universal_subbundle(Gr(2, 5));

julia> degree(cohomology(U)[0])
0

julia> degree(cohomology(dual(U))[0])
5
```
"""
function universal_subbundle(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal subbundle")
  )
  return dual(CompletelyReducibleBundle(X, fundamental_weight(dynkin_type(X), 1)))
end

"""
    universal_quotient_bundle(X::PartialFlagVariety) -> Union{CompletelyReducibleBundle, FilteredBundle}

The universal quotient bundle ``\\mathcal{Q}`` on a generalized Grassmannian.

On ``\\mathrm{Gr}(k, n) = \\mathrm{A}_{n-1}/P_k`` (type ``\\mathrm{A}``),
``\\mathcal{Q}`` is the completely reducible equivariant bundle with highest
weight ``\\omega_{n-1}`` — geometrically the rank-``(n-k)`` quotient
``\\mathbb{C}^n / \\mathcal{U}``.

For isotropic Grassmannians (types ``\\mathrm{B}``, ``\\mathrm{C}``,
``\\mathrm{D}``), there is no literal quotient of the standard representation
by ``\\mathcal{U}``; instead ``\\mathcal{Q}`` is taken to be
``(\\mathcal{U}^\\perp)^\\vee``, dual to the orthogonal complement of
``\\mathcal{U}``. When the residual bundle ``\\mathcal{R}`` vanishes
(Lagrangian / spinor cases) ``\\mathcal{Q} \\cong \\mathcal{U}^\\vee``;
otherwise ``\\mathcal{Q}`` is a filtered bundle with graded pieces
``\\mathcal{R}^\\vee`` and ``\\mathcal{U}^\\vee``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> Q = universal_quotient_bundle(Gr(2, 5));

julia> rank_bundle(Q)
3

julia> degree(cohomology(Q)[0])  # H⁰(Q) is the standard representation of A_4
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

``\\mathcal{R}`` vanishes on Lagrangian Grassmannians (the spinor / maximal isotropic
cases), and is acyclic on all isotropic Grassmannians since its highest weights
are not dominant for the parabolic:
```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(residual_bundle(SGr(3, 6)))  # Lagrangian: R vanishes
0

julia> degree(cohomology(residual_bundle(OGr(2, 7)))[0])  # acyclic on OGr(2, 7)
0
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

```jldoctest
julia> using PartialFlagVarieties

julia> X = quadric(5);  # B_3/P_1

julia> degree(cohomology(spinor_bundle(X))[0])  # spin representation of Spin(7)
8
```
"""
function spinor_bundle(X::PartialFlagVariety)
  _check_spinor_domain(X)
  DT = dynkin_type(X)
  R = rank(DT)
  # OGr(k, 2n+1): single spinor bundle at ω_n.
  DT <: TypeB && return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  # DT <: TypeD; OGr(k, 2n) and D_n/P_{n-1, n} carry both half-spinors as a direct sum.
  return CompletelyReducibleBundle(
    X, [fundamental_weight(DT, R - 1), fundamental_weight(DT, R)]
  )
end

function spinor_bundle(X::PartialFlagVariety, half::Symbol)
  _check_spinor_domain(X)
  DT = dynkin_type(X)
  R = rank(DT)

  if DT <: TypeB
    half in (:plus, :minus) &&
      @warn "Type B has a single spinor bundle; ignoring half=$half"
    return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  end

  half === :plus && return CompletelyReducibleBundle(X, fundamental_weight(DT, R - 1))
  half === :minus && return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
  throw(ArgumentError("half must be :plus or :minus, got :$half"))
end

# Spinor bundles are defined on B_n/P_k, D_n/P_k, and the two-marked D_n/P_{n-1, n}.
# The two-marked case is the (n-1)-isotropic Grassmannian (Picard rank 2) — see
# Frassineti–Manivel, arXiv:2605.28712, §1.
function _check_spinor_domain(X::PartialFlagVariety)
  is_orthogonal_grassmannian(X) && return nothing
  DT = dynkin_type(X)
  if DT <: TypeD
    marked = marked_nodes(X)
    R = rank(DT)
    length(marked) == 2 && Set(marked) == Set((R - 1, R)) && return nothing
  end
  throw(
    ArgumentError(
      "spinor_bundle requires B_n/P_k, D_n/P_k, or D_n/P_{n-1, n}"
    ),
  )
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Partial flag varieties Fl(d₁,...,dₖ; n)
# ═══════════════════════════════════════════════════════════════════════════════
# TODO: Extend to isotropic types. For a generalized Grassmannian of type B, C, or D
# the building blocks are universal_subbundle(X) and residual_bundle(X). For multi-step
# isotropic flags, compute the graded pieces U_i/U_{i-1} via weight differences
# analogously to the type A case below.

"""
    universal_subbundle(X::PartialFlagVariety, i::Int) -> Bundle

The ``i``-th universal subbundle ``\\mathcal{U}_i`` on a partial flag variety
``\\mathrm{Fl}(d_1, \\ldots, d_k; n)``, selected from the filtration
``0 \\subset \\mathcal{U}_1 \\subset \\cdots \\subset \\mathcal{U}_k``.

`i == 1` returns [`universal_subbundle(X)`](@ref). For `i > 1`, equivalent to
`universal_subbundles(X)[i]`.

Currently `i > 1` is only supported for type ``\\mathrm{A}`` partial flag
varieties.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, [1, 2]);  # Fl(1,2; 4)

julia> rank_bundle(universal_subbundle(X, 1))
1

julia> rank_bundle(universal_subbundle(X, 2))
2
```
"""
function universal_subbundle(X::PartialFlagVariety, i::Int)
  i == 1 && return universal_subbundle(X)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal subbundle")
  )
  # TODO: For isotropic Grassmannians (types B, C, D) with i == 2, return U^⊥ as
  # dual(universal_quotient_bundle(X)), reflecting 0 → U → U^⊥ → Res → 0.
  # When Res is zero (Lagrangian/spinor cases), U^⊥ = U, so dual(Q) recovers i == 1.
  k = length(marked_nodes(X))
  (1 ≤ i ≤ k) ||
    throw(ArgumentError("i must be between 1 and length(marked_nodes(X)) = $k"))
  return universal_subbundles(X)[i]
end

"""
    tautological_bundles(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

The graded pieces ``\\mathcal{U}_1, \\mathcal{U}_2 / \\mathcal{U}_1, \\ldots,
\\mathcal{U}_k / \\mathcal{U}_{k-1}`` of the tautological filtration on a
partial flag variety ``\\mathrm{Fl}(d_1, \\ldots, d_k; n)``, returned as a
vector of completely reducible bundles.

These are the convenient building blocks; for the filtration steps themselves
(as `FilteredBundle`s), use [`universal_subbundles`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, [1, 2]);  # Fl(1,2; 4)

julia> U = tautological_bundles(X);

julia> rank_bundle.(U)
2-element Vector{Int64}:
 1
 1
```
"""
function tautological_bundles(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have well-defined tautological bundles")
  )
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  # TODO: For isotropic types B, C, D, the graded pieces of the tautological filtration are
  # universal_subbundle(X) and residual_bundle(X) for Grassmannians, and more generally
  # the U_i/U_{i-1} computed from the weight differences at each marked node.
  DT <: TypeA || throw(
    ArgumentError(
      "tautological_bundles not implemented for non-type A partial flag varieties"
    ),
  )
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
end

"""
    universal_subbundles(X::PartialFlagVariety) -> Vector{Bundle}

The full flag of universal subbundles on a partial flag variety
``\\mathrm{Fl}(d_1, \\ldots, d_k; n)``:
``0 \\subset \\mathcal{U}_1 \\subset \\cdots \\subset \\mathcal{U}_k``,
returned as a vector. The first element is the rank-``d_1`` subbundle
``\\mathcal{U}_1``; the ``i``-th element (for ``i > 1``) is the `FilteredBundle`
with graded pieces ``\\mathcal{U}_1, \\mathcal{U}_2 / \\mathcal{U}_1, \\ldots,
\\mathcal{U}_i / \\mathcal{U}_{i-1}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, [1, 2]);  # Fl(1,2; 4)

julia> subs = universal_subbundles(X);

julia> rank_bundle.(subs)
2-element Vector{Int64}:
 1
 2
```
"""
function universal_subbundles(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have well-defined universal subbundles")
  )
  DT = dynkin_type(X)
  # TODO: For isotropic Grassmannians (types B, C, D), the filtration steps are
  # [U] and [U^⊥ = FilteredBundle(X, [universal_subbundle(X), residual_bundle(X)])].
  # For isotropic partial flags, extend similarly to the type A implementation above.
  DT <: TypeA || throw(
    ArgumentError(
      "universal_subbundles currently only implemented for type A partial flag varieties"
    ),
  )
  marked = marked_nodes(X)
  bundles = Vector{Bundle}()
  U = tautological_bundles(X)
  push!(bundles, U[1])
  for i in 2:length(marked)
    push!(bundles, FilteredBundle(X, U[1:i]))
  end
  return bundles
end
