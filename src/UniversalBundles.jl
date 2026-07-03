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
export S, Q  # shorthands for universal_subbundle, universal_quotient_bundle

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

``\\mathcal{U}`` has no dominant weights (its ``\\rho``-shifted weights are
singular), so ``\\mathrm{H}^0(\\mathcal{U}) = 0``,
while ``\\mathrm{H}^0(\\mathcal{U}^\\vee)`` recovers the standard representation of the
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
    S(X::PartialFlagVariety) -> CompletelyReducibleBundle
    S(X::PartialFlagVariety, i::Int) -> Bundle

Shorthand for [`universal_subbundle`](@ref): the tautological subbundle
``\\mathcal{S}`` on `X`. On a multi-step flag, `S(X, i)` is the `i`-th universal
subbundle.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(S(Gr(2, 5)))
2
```
"""
S(X::PartialFlagVariety) = universal_subbundle(X)
S(X::PartialFlagVariety, i::Int) = universal_subbundle(X, i)

"""
    Q(X::PartialFlagVariety) -> Union{CompletelyReducibleBundle, FilteredBundle}

Shorthand for [`universal_quotient_bundle`](@ref): the universal quotient bundle
``\\mathcal{Q}`` on a generalized Grassmannian `X`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(Q(Gr(2, 5)))
3
```
"""
Q(X::PartialFlagVariety) = universal_quotient_bundle(X)

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
cases); on the isotropic Grassmannians below it is acyclic because the
``\\rho``-shifts of its weights are singular (by Bott, non-dominance alone
would only move cohomology to higher degree):
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

The spinor bundle on an orthogonal Grassmannian ``\\mathrm{OGr}(k, m)``. The two
methods are mutually exclusive by Dynkin type.

For type ``\\mathrm{B}_n/P_k`` (``\\mathrm{OGr}(k, 2n+1)``) there is a single
spinor bundle ``\\mathcal{S}`` corresponding to ``\\omega_n``; use the
one-argument `spinor_bundle(X)`.

For type ``\\mathrm{D}_n/P_k`` (``\\mathrm{OGr}(k, 2n)``) there are two
half-spinor bundles ``\\mathcal{S}^+`` and ``\\mathcal{S}^-`` corresponding to
``\\omega_{n-1}`` and ``\\omega_n`` respectively; select one with
`spinor_bundle(X, :plus)` or `spinor_bundle(X, :minus)`.

Calling `spinor_bundle(X)` on type ``\\mathrm{D}``, or `spinor_bundle(X, half)`
on type ``\\mathrm{B}``, throws an `ArgumentError`.

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
  DT <: TypeB || throw(
    ArgumentError(
      "spinor_bundle(X) is only defined for type B (a single spinor bundle); " *
      "type D has two half-spinors — use spinor_bundle(X, :plus) or spinor_bundle(X, :minus)",
    ),
  )
  # OGr(k, 2n+1): single spinor bundle at ω_n.
  return CompletelyReducibleBundle(X, fundamental_weight(DT, rank(DT)))
end

function spinor_bundle(X::PartialFlagVariety, half::Symbol)
  _check_spinor_domain(X)
  DT = dynkin_type(X)
  DT <: TypeD || throw(
    ArgumentError(
      "spinor_bundle(X, half) is only defined for type D (two half-spinors); " *
      "type B has a single spinor bundle — use spinor_bundle(X)",
    ),
  )
  half in (:plus, :minus) ||
    throw(ArgumentError("half must be :plus or :minus, got :$half"))
  R = rank(DT)
  # OGr(k, 2n): ω_{n-1} (plus) and ω_n (minus).
  half === :plus && return CompletelyReducibleBundle(X, fundamental_weight(DT, R - 1))
  return CompletelyReducibleBundle(X, fundamental_weight(DT, R))
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

"""
    universal_subbundle(X::PartialFlagVariety, i::Int) -> Bundle

The ``i``-th universal subbundle on a partial flag variety, selected from
[`universal_subbundles`](@ref). `i == 1` returns
[`universal_subbundle(X)`](@ref).

On an isotropic generalized Grassmannian (one marked node, type
``\\mathrm{B}``, ``\\mathrm{C}``, or ``\\mathrm{D}``), `i == 2` returns the
orthogonal complement ``\\mathcal{U}^\\perp`` — there is a canonical second
subbundle even though only one node is marked.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, [1, 2]);  # Fl(1,2; 4)

julia> rank_bundle(universal_subbundle(X, 1))
1

julia> rank_bundle(universal_subbundle(X, 2))
2
```

```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle(universal_subbundle(OGr(3, 9), 2))  # U^⊥, rank n-k = 6
6
```
"""
function universal_subbundle(X::PartialFlagVariety, i::Int)
  i == 1 && return universal_subbundle(X)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have a well-defined universal subbundle")
  )
  bundles = universal_subbundles(X)
  (1 <= i <= length(bundles)) || throw(
    ArgumentError(
      "i must be between 1 and length(universal_subbundles(X)) = $(length(bundles))"
    ),
  )
  return bundles[i]
end

"""
    tautological_bundles(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

The graded pieces of the tautological filtration on a partial flag variety.

For ``\\mathrm{Fl}(d_1, \\ldots, d_k; n)`` (type ``\\mathrm{A}``) and
multi-step isotropic partial flags (type ``\\mathrm{B}``, ``\\mathrm{C}``,
``\\mathrm{D}``), returns ``k`` pieces
``\\mathcal{U}_1, \\mathcal{U}_2/\\mathcal{U}_1, \\ldots, \\mathcal{U}_k/\\mathcal{U}_{k-1}``.

For an isotropic **generalized Grassmannian** (one marked node, type
``\\mathrm{B}``, ``\\mathrm{C}``, or ``\\mathrm{D}``), returns the two graded
pieces ``[\\mathcal{U}, \\mathcal{R}]`` of the natural filtration
``0 \\subset \\mathcal{U} \\subset \\mathcal{U}^\\perp``. In Lagrangian /
spinor cases ``\\mathcal{R}`` is the zero bundle.

The two-marked spinorial variety ``\\mathrm{D}_n/P_{n-1, n}`` is also
supported: it returns the two graded pieces of the unique
``(n-1) \\subset n`` isotropic flag (ranks ``[n-1, 1]``).

# References

- Pragacz–Ratajski, "Formulas for Lagrangian and orthogonal degeneracy loci",
  *Compositio Math.* **107** (1997), 11–87 — tautological bundles on
  isotropic flag varieties.
- Buch–Kresch–Tamvakis, "Quantum Pieri rules for isotropic Grassmannians",
  *Invent. Math.* **178** (2009), 345–405.
- Bourbaki, *Groupes et algèbres de Lie*, Ch. VI, Planche IV — fundamental
  weights of ``\\mathrm{D}_n``, used in the spinor-boundary correction below.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, [1, 2]);  # Fl(1, 2; 4)

julia> rank_bundle.(tautological_bundles(X))
2-element Vector{Int64}:
 1
 1
```

```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle.(tautological_bundles(OGr(3, 9)))  # [U, R] on OGr(3, 9)
2-element Vector{Int64}:
 3
 3
```
"""
function tautological_bundles(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have well-defined tautological bundles")
  )
  DT = dynkin_type(X)
  marked = marked_nodes(X)

  if DT <: TypeA
    return _tautological_pieces_typeA(X, DT, marked)
  end
  DT <: Union{TypeB,TypeC,TypeD} || throw(
    ArgumentError("tautological_bundles not implemented for $(DT) partial flag varieties")
  )

  # Isotropic generalized Grassmannian: the natural filtration is 0 ⊂ U ⊂ U^⊥,
  # so the graded pieces are [U, R] (Pragacz–Ratajski 1997, §1).
  is_generalized_grassmannian(X) &&
    return [universal_subbundle(X), residual_bundle(X)]

  # Multi-step isotropic partial flag: same Levi-rep highest-weight pattern as
  # type A, with one D-specific correction at the spinor boundary.
  return _tautological_pieces_isotropic(X, DT, marked)
end

function _tautological_pieces_typeA(X, DT, marked)
  bundles = Vector{CompletelyReducibleBundle}()
  push!(bundles, dual(CompletelyReducibleBundle(X, fundamental_weight(DT, 1))))
  for i in 2:length(marked)
    push!(
      bundles,
      CompletelyReducibleBundle(
        X,
        -fundamental_weight(DT, marked[i]) + fundamental_weight(DT, marked[i] - 1),
      ),
    )
  end
  return bundles
end

function _tautological_pieces_isotropic(X, DT, marked)
  R = rank(DT)
  bundles = Vector{CompletelyReducibleBundle}()
  push!(bundles, dual(CompletelyReducibleBundle(X, fundamental_weight(DT, 1))))
  for i in 2:length(marked)
    m = marked[i]
    # Highest weight of U_{m_i} / U_{m_{i-1}} as a Levi rep. Away from the
    # type-D spinor boundary, the type-A formula ω_{m-1} - ω_m = L_m
    # (Bourbaki, Groupes et algèbres de Lie, Ch. VI, Plate IV) gives the
    # standard rep of the i-th GL Levi factor.
    #
    # In type D at m = n - 1, ω_{n-1} = (L_1 + ⋯ + L_{n-1} - L_n)/2 and
    # ω_n = (L_1 + ⋯ + L_n)/2 are half-sums, so ω_{n-2} - ω_{n-1} is
    # half-integer and picks up a spinor Levi character with the wrong rank.
    # The correction ω_n - ω_{n-1} = L_n is the genuine basis-vector weight,
    # giving rank n - m_{i-1} — matching U_{n-1} = max isotropic of dim n in
    # the "minus" spinor family (Manivel, SIGMA 5 (2009) 078). At m = n the
    # type-A formula already produces ω_{n-1} - ω_n = -L_n.
    ω = if DT <: TypeD && m == R - 1
      fundamental_weight(DT, R) - fundamental_weight(DT, R - 1)
    else
      fundamental_weight(DT, m - 1) - fundamental_weight(DT, m)
    end
    push!(bundles, CompletelyReducibleBundle(X, ω))
  end
  return bundles
end

"""
    universal_subbundles(X::PartialFlagVariety) -> Vector{Bundle}

The flag of universal subbundles on a partial flag variety, returned as a
vector. The first element is ``\\mathcal{U}_1`` as a
`CompletelyReducibleBundle`; later elements are `FilteredBundle`s assembled
from the graded pieces returned by [`tautological_bundles`](@ref).

### Type A
For ``\\mathrm{Fl}(d_1, \\ldots, d_k; n)`` the returned vector is the
nested flag of tautological subbundles
``\\mathcal{U}_1 \\subset \\mathcal{U}_2 \\subset \\cdots \\subset
\\mathcal{U}_k``, with ``\\mathrm{rank}(\\mathcal{U}_i) = d_i``.

### Isotropic generalized Grassmannian (one marked node, B/C/D)
The natural filtration is ``0 \\subset \\mathcal{U} \\subset
\\mathcal{U}^\\perp``, so the function returns
``[\\mathcal{U}, \\mathcal{U}^\\perp]``. In Lagrangian and spinor cases
``\\mathcal{R} = 0`` and the two entries coincide.

### Orthogonal / symplectic multi-step partial flag (k ≥ 2 marked nodes)
The marked nodes ``m_1 < m_2 < \\cdots < m_k`` parameterise nested
isotropic subspaces; entry ``i`` is the corresponding tautological
subbundle ``\\mathcal{U}_{m_i}``. In types B and C and in type D away
from the spinor nodes, ``\\mathrm{rank}(\\mathcal{U}_{m_i}) = m_i``.

In type D, the parabolics ``P_{n-1}`` and ``P_n`` correspond to the two
families of *maximal* isotropic subspaces (both of dimension ``n``). When a
last marked node ``m_k \\in \\{n-1, n\\}`` appears in a multi-step flag,
the corresponding ``\\mathcal{U}_{m_k}`` is the rank-``n`` maximal
isotropic from that spinor family — *not* a rank-``(n-1)`` subspace. This
follows the convention of the underlying generalized Grassmannian: e.g.
`OGr(3, 8) = D_4 / P_3` has `rank_bundle(universal_subbundle(X)) == 4`.

The two-marked spinorial variety ``\\mathrm{D}_n / P_{n-1, n}`` (the
``(n-1)``-isotropic Grassmannian, Picard rank 2 — see Frassineti–Manivel,
arXiv:2605.28712, §1) is supported: it returns
``[\\mathcal{U}_{n-1}, \\mathcal{U}_n]`` with ranks ``[n-1, n]``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = flag_variety(4, [1, 2]);  # Fl(1, 2; 4)

julia> rank_bundle.(universal_subbundles(X))
2-element Vector{Int64}:
 1
 2
```

```jldoctest
julia> using PartialFlagVarieties

julia> rank_bundle.(universal_subbundles(OGr(3, 9)))  # [U, U^⊥]
2-element Vector{Int64}:
 3
 6
```

```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeD{4}, (1, 3));  # last marked at the spinor node

julia> rank_bundle.(universal_subbundles(X))
2-element Vector{Int64}:
 1
 4
```
"""
function universal_subbundles(X::PartialFlagVariety)
  is_exceptional_type(X) && throw(
    ArgumentError("exceptional types do not have well-defined universal subbundles")
  )
  DT = dynkin_type(X)
  DT <: Union{TypeA,TypeB,TypeC,TypeD} || throw(
    ArgumentError("universal_subbundles not implemented for $(DT) partial flag varieties")
  )
  bundles = Vector{Bundle}()
  U = tautological_bundles(X)
  push!(bundles, U[1])
  for i in 2:length(U)
    push!(bundles, FilteredBundle(X, U[1:i]))
  end
  return bundles
end
