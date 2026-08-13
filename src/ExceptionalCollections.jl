# ═══════════════════════════════════════════════════════════════════════════════
#  ExceptionalCollections — exceptional objects and collections on G/P
#
#  Implements predicates for exceptionality and standard constructions:
#
#   Predicates
#   - is_exceptional            check Ext*(E, E) = k in degree 0
#   - is_exceptional_pair       check RHom(F, E) = 0
#   - is_strong_exceptional_pair  check Hom(E,F) free + Ext^{>0}(E,F)=0
#   - is_exceptional_sequence   pairwise exceptionality
#   - is_strong_exceptional_sequence  pairwise strong exceptionality
#   - is_full_exceptional_sequence    check length = χ(X)
#
#   Standard collections
#   - beilinson_collection      O, O(1), ..., O(n) on ℙⁿ
#   - kapranov_collection       Kapranov's collection on Q^n and Gr(k,n)
#   - schur_functor             Σ^α(E) for an irreducible E (Type A, Gr)
#   - kapranov_bundles_grassmannian  all Σ^α U^∨ bundles on Gr(k,n)
# ═══════════════════════════════════════════════════════════════════════════════

export is_exceptional, is_exceptional_pair, is_strong_exceptional_pair
export is_exceptional_sequence, is_strong_exceptional_sequence
export is_full_exceptional_sequence
export beilinson_collection, beilinson_collection_dual
export kapranov_collection
export schur_functor, kapranov_bundles_grassmannian

const _SEQUENCE_STATUS_CACHE = let b = _default_cache_budget()
  LRU{Tuple{Vector{CompletelyReducibleBundle},CompletelyReducibleBundle},Tuple{Bool,Bool}}(;
    maxsize=_cache_maxsize(b, _DEFAULT_STRUCTURAL_FRAC * 0.1),
    by=Base.summarysize,
  )
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Exceptionality predicates
# ═══════════════════════════════════════════════════════════════════════════════

"""
    is_exceptional(E::CompletelyReducibleBundle) -> Bool

Check whether `E` is an exceptional object in `D^b(G/P)`.

An object `E` is **exceptional** if:
- ``\\operatorname{Hom}(\\mathcal{E}, \\mathcal{E}) = \\mathbb{C}``  (i.e., ``\\mathrm{h}^0(\\mathcal{E}^\\vee \\otimes \\mathcal{E}) = 1``), and
- ``\\operatorname{Ext}^i(\\mathcal{E}, \\mathcal{E}) = 0`` for all ``i > 0``  (i.e., ``\\mathrm{h}^i(\\mathcal{E}^\\vee \\otimes \\mathcal{E}) = 0``).

Uses Borel–Weil–Bott to compute cohomology of ``\\mathcal{E}^\\vee \\otimes \\mathcal{E}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> is_exceptional(structure_sheaf(X))
true

julia> is_exceptional(line_bundle(X, 2))
true
```
"""
function is_exceptional(E::CompletelyReducibleBundle)
  H = cohomology(dual(E) ⊗ E)
  H[0] == 1 && all(H[i] == 0 for i in 1:(H.max_degree))
end

"""
    is_exceptional_pair(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle) -> Bool

Check whether `(E, F)` is an **exceptional pair**: ``\\operatorname{RHom}(\\mathcal{F}, \\mathcal{E}) = 0``,
i.e., ``\\mathrm{H}^i(\\mathrm{G}/\\mathrm{P},\\, \\mathcal{F}^\\vee \\otimes \\mathcal{E}) = 0`` for all ``i \\geq 0``.

Note: this is the orthogonality condition; the pair is ordered so that
``(\\mathcal{E}_i, \\mathcal{E}_j)`` is an exceptional pair for ``i < j`` in an exceptional
sequence ``\\langle \\mathcal{E}_1, \\ldots, \\mathcal{E}_n \\rangle``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> O0 = structure_sheaf(X);

julia> O1 = line_bundle(X, 1);

julia> is_exceptional_pair(O0, O1)
true
```
"""
function is_exceptional_pair(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  # Ext^i(F, E) = H^i(F^∨ ⊗ E)
  FvE = dual(F) ⊗ E
  return iszero(cohomology(FvE))
end

"""
    is_strong_exceptional_pair(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle) -> Bool

Check whether `(E, F)` is a **strong exceptional pair**:
- ``\\operatorname{Ext}^i(\\mathcal{F}, \\mathcal{E}) = 0`` for all ``i`` (i.e., ``(\\mathcal{E}, \\mathcal{F})`` is an exceptional pair), and
- ``\\operatorname{Ext}^i(\\mathcal{E}, \\mathcal{F}) = 0`` for all ``i > 0`` (only ``\\operatorname{Hom}(\\mathcal{E}, \\mathcal{F})`` can be nonzero).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> O0 = structure_sheaf(X);

julia> O1 = line_bundle(X, 1);

julia> is_strong_exceptional_pair(O0, O1)
true
```
"""
function is_strong_exceptional_pair(
  E::CompletelyReducibleBundle, F::CompletelyReducibleBundle
)
  # Check RHom(F, E) = 0
  is_exceptional_pair(E, F) || return false
  # Check Ext^{>0}(E, F) = 0, i.e., H^i(E^∨ ⊗ F) = 0 for i > 0
  _has_no_positive_exts(E, F)
end

function _has_no_positive_exts(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  H = cohomology(dual(E) ⊗ F)
  all(H[i] == 0 for i in 1:(H.max_degree))
end

"""
    is_exceptional_sequence(Es::Vector{<:CompletelyReducibleBundle}) -> Bool

Check whether `Es = [E₁, ..., Eₙ]` is an **exceptional sequence**:
- Each ``\\mathcal{E}_i`` is exceptional (``\\operatorname{Ext}^\\vee(\\mathcal{E}_i, \\mathcal{E}_i) = \\mathbb{C}``), and
- ``(\\mathcal{E}_i, \\mathcal{E}_j)`` is an exceptional pair (``\\operatorname{RHom}(\\mathcal{E}_j, \\mathcal{E}_i) = 0``) for all ``i < j``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> Es = [line_bundle(X, k) for k in 0:3];

julia> is_exceptional_sequence(Es)
true
```
"""
function is_exceptional_sequence(Es::Vector{<:CompletelyReducibleBundle})
  n = length(Es)
  all(is_exceptional, Es) &&
    all(is_exceptional_pair(Es[i], Es[j]) for i in 1:n for j in (i + 1):n)
end

"""
    is_strong_exceptional_sequence(Es::Vector{<:CompletelyReducibleBundle}) -> Bool

Check whether `Es = [E₁, ..., Eₙ]` is a **strong exceptional sequence**:
- Each ``\\mathcal{E}_i`` is exceptional, and
- ``(\\mathcal{E}_i, \\mathcal{E}_j)`` is a strong exceptional pair for all ``i < j``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> Es = [line_bundle(X, k) for k in 0:3];

julia> is_strong_exceptional_sequence(Es)
true
```
"""
function is_strong_exceptional_sequence(Es::Vector{<:CompletelyReducibleBundle})
  n = length(Es)
  all(is_exceptional, Es) &&
    all(is_strong_exceptional_pair(Es[i], Es[j]) for i in 1:n for j in (i + 1):n)
end

"""
    is_full_exceptional_sequence(Es::Vector{<:CompletelyReducibleBundle},
                                  X::PartialFlagVariety) -> Bool

Check whether `Es` is a **full exceptional sequence** on `X`.

Necessary condition used: the sequence is exceptional and its length equals
the Euler characteristic ``\\chi(X) = \\sum_i b_{2i}`` (which for ``\\mathrm{G}/\\mathrm{P}`` equals
the total number of Schubert cells / the number of T-fixed points).

For ``\\mathrm{G}/\\mathrm{P}``, which has only even cohomology, ``\\chi(X) = \\sum b_i``.

!!! note
    This is a necessary, not sufficient, condition for fullness.
    It verifies the K-theoretic Grothendieck group rank (a necessary
    numerical condition), not fullness in the derived category.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> Es = [line_bundle(X, k) for k in 0:3];

julia> is_full_exceptional_sequence(Es, X)
true
```
"""
function is_full_exceptional_sequence(
  Es::Vector{<:CompletelyReducibleBundle},
  X::PartialFlagVariety,
)
  χ = euler_characteristic(X)
  length(Es) == χ || return false
  return is_exceptional_sequence(Es)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Beilinson collection on ℙⁿ
# ═══════════════════════════════════════════════════════════════════════════════

"""
    beilinson_collection(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

The **Beilinson exceptional collection** on projective space ``\\mathbb{P}^n``:

``\\langle \\mathcal{O}, \\mathcal{O}(1), \\ldots, \\mathcal{O}(n) \\rangle``

This is a full strong exceptional collection on ``\\mathbb{P}^n``.

`X` must be a **type-A** projective space, i.e. ``\\mathrm{A}_n/\\mathrm{P}_1`` or
``\\mathrm{A}_n/\\mathrm{P}_n``. Beilinson's construction is intrinsically type A, so the
other models of ``\\mathbb{P}^n`` recognised by `is_projective_space`
(``\\mathrm{C}_n/\\mathrm{P}_1``, ``\\mathrm{B}_2/\\mathrm{P}_2``, ``\\mathrm{D}_3/\\mathrm{P}_2``,
``\\mathrm{D}_3/\\mathrm{P}_3``) are rejected.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> Es = beilinson_collection(X);

julia> length(Es)
4

julia> all(is_exceptional, Es)
true
```
"""
function beilinson_collection(X::PartialFlagVariety)
  is_projective_space(X) || throw(ArgumentError(
    "beilinson_collection requires a projective space ℙⁿ"))
  dynkin_type(X) <: TypeA || throw(
    ArgumentError(
      "beilinson_collection is only implemented for type-A ℙⁿ; got a non-type-A " *
      "model of projective space (Cₙ/P₁, B₂/P₂, D₃/P₂, or D₃/P₃)"),
  )
  n = dimension(X)   # = rank(DT) = n for ℙⁿ
  [line_bundle(X, k) for k in 0:n]
end

"""
    beilinson_collection_dual(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

The **dual Beilinson exceptional collection** on ``\\mathbb{P}^n``:

``\\langle \\Omega^n(n), \\Omega^{n-1}(n-1), \\ldots, \\Omega^1(1), \\mathcal{O} \\rangle``

where ``\\Omega^k(k) = \\bigwedge\\nolimits^k \\Omega^1_{\\mathbb{P}^n} \\otimes \\mathcal{O}(k)``.

`X` must be a **type-A** projective space, i.e. ``\\mathrm{A}_n/\\mathrm{P}_1`` or
``\\mathrm{A}_n/\\mathrm{P}_n``. The construction is intrinsically type A: for the other
models of ``\\mathbb{P}^n`` recognised by `is_projective_space`
(``\\mathrm{C}_n/\\mathrm{P}_1``, ``\\mathrm{B}_2/\\mathrm{P}_2``, ``\\mathrm{D}_3/\\mathrm{P}_2``,
``\\mathrm{D}_3/\\mathrm{P}_3``) the homogeneous cotangent bundle does not reproduce
``\\Omega^k(k)``, so they are rejected.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> Es = beilinson_collection_dual(X);

julia> length(Es)
4

julia> all(is_exceptional, Es)
true
```
"""
function beilinson_collection_dual(X::PartialFlagVariety)
  is_projective_space(X) || throw(ArgumentError(
    "beilinson_collection_dual requires a projective space ℙⁿ"))
  dynkin_type(X) <: TypeA || throw(
    ArgumentError(
      "beilinson_collection_dual is only implemented for type-A ℙⁿ; got a " *
      "non-type-A model of projective space (Cₙ/P₁, B₂/P₂, D₃/P₂, or D₃/P₃)"),
  )
  n = dimension(X)
  Omega = cotangent_bundle(X)
  # Omega^k(k) = ∧^k Omega ⊗ O(k) for k = n, …, 1, followed by O.
  collection = CompletelyReducibleBundle[
    twist(exterior_power(Omega, k), 1, k) for k in n:-1:1
  ]
  push!(collection, structure_sheaf(X))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Kapranov collection on quadrics
# ═══════════════════════════════════════════════════════════════════════════════

"""
    kapranov_collection(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

The **Kapranov exceptional collection** on a smooth quadric ``Q^n``.

- **Odd quadric** ``Q^{2m-1}`` (``\\mathrm{B}_m/\\mathrm{P}_1``):
  ``\\langle \\mathcal{O},\\; \\Sigma,\\; \\mathcal{O}(1),\\; \\ldots,\\; \\mathcal{O}(n-1) \\rangle``
  where ``\\Sigma`` is the spinor bundle.
- **Even quadric** ``Q^{n}`` (``\\mathrm{D}_{n/2+1}/\\mathrm{P}_1``):
  ``\\langle \\mathcal{O},\\; \\Sigma^+,\\; \\Sigma^-,\\; \\mathcal{O}(1),\\; \\ldots,\\; \\mathcal{O}(n-1) \\rangle``
  where ``\\Sigma^\\pm`` are the two half-spinor bundles.

`X` must be a quadric (``\\mathrm{B}_m/\\mathrm{P}_1`` or ``\\mathrm{D}_m/\\mathrm{P}_1``).

# References
- M. Kapranov, *On the derived categories of coherent sheaves on some homogeneous
  spaces*, Invent. Math. 92 (1988), 479–508.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = quadric(5);  # Q^5 = B_3/P_1

julia> Es = kapranov_collection(X);

julia> length(Es)
6

julia> is_full_exceptional_sequence(Es, X)
true
```
"""
function kapranov_collection(X::PartialFlagVariety)
  is_quadric(X) || throw(ArgumentError(
    "kapranov_collection requires a quadric (B_m/P_1 or D_m/P_1)"
  ))
  DT = dynkin_type(X)
  n = dimension(X)
  result = CompletelyReducibleBundle[]

  if DT <: TypeB
    # Q^{2m-1}: ⟨O, Σ, O(1), ..., O(n-1)⟩  (n+1 bundles = χ(Q^n))
    push!(result, structure_sheaf(X))
    push!(result, spinor_bundle(X))
    for k in 1:(n - 1)
      push!(result, line_bundle(X, k))
    end
  elseif DT <: TypeD
    # Q^n (n even): ⟨O, Σ⁺, Σ⁻, O(1), ..., O(n-1)⟩  (n+2 bundles = χ(Q^n))
    push!(result, structure_sheaf(X))
    push!(result, spinor_bundle(X, :plus))
    push!(result, spinor_bundle(X, :minus))
    for k in 1:(n - 1)
      push!(result, line_bundle(X, k))
    end
  end
  result
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Schur functor on Type-A Grassmannians
# ═══════════════════════════════════════════════════════════════════════════════

"""
    schur_functor(X::PartialFlagVariety, partition::Vector{<:Integer})
        -> CompletelyReducibleBundle

The Schur functor ``\\Sigma^\\alpha \\mathcal{U}^\\vee`` on the Grassmannian
``\\mathrm{Gr}(k, n) = \\mathrm{A}_{n-1}/\\mathrm{P}_k``, where ``\\alpha`` is a Young diagram
with rows given by `partition`.

The partition must fit in a ``k \\times (n-k)`` box:
rows ``\\alpha_1 \\geq \\alpha_2 \\geq \\ldots \\geq \\alpha_k \\geq 0``, all ``\\leq n-k``.

The ambient weight in ``\\mathrm{A}_{n-1}`` is
```math
\\lambda = \\sum_{i=1}^{k-1} (\\alpha_i - \\alpha_{i+1})\\,\\omega_i + \\alpha_k\\,\\omega_k
```
where ``\\omega_1, \\ldots, \\omega_k`` are the first ``k`` fundamental weights of
``\\mathrm{A}_{n-1}`` (corresponding to the semisimple ``\\mathrm{A}_{k-1}``-part of the Levi at nodes
``1, \\ldots, k-1`` and the central character at the marked node ``k``).

`X` must be a Type-A generalized Grassmannian.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> # Σ^{(1,0)} U^∨ = U^∨  (rank 2)
       E = schur_functor(X, [1, 0]);

julia> rank(E)
2

julia> # Σ^{(1,1)} U^∨ = det(U^∨) = O(1)  (rank 1)
       L1 = schur_functor(X, [1, 1]);

julia> rank(L1)
1
```
"""
function schur_functor(
  X::PartialFlagVariety,
  partition::Vector{<:Integer},
)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  DT <: TypeA || throw(ArgumentError(
    "schur_functor is only implemented for Type A Grassmannians"
  ))
  length(marked) == 1 || throw(
    ArgumentError(
      "schur_functor requires a generalized Grassmannian (one marked node)"
    ),
  )

  k = Int(marked[1])     # Gr(k, n) → marked node k
  n = rank(DT) + 1       # = n  (since rank(A_{n-1}) = n-1)
  nk = n - k             # n - k

  # Pad/truncate partition to length k
  α = zeros(Int, k)
  for i in 1:min(length(partition), k)
    α[i] = Int(partition[i])
  end

  # Validate
  issorted(α; rev=true) || throw(ArgumentError(
    "Partition must be weakly decreasing, got $partition"
  ))
  all(>=(0), α) || throw(ArgumentError("Partition entries must be ≥ 0"))
  α[1] <= nk || throw(
    ArgumentError(
      "Partition $(partition) does not fit in a $k × $nk box (max row length = $(α[1]) > $nk)"
    ),
  )

  # Compute ambient weight using the GL(k) Schur functor formula.
  # The Levi of P_k in A_{n-1} has semisimple part A_{k-1} sitting at nodes 1,...,k-1,
  # with central character at the marked node k.
  # For partition α = (α₁ ≥ ... ≥ α_k ≥ 0), the ambient weight is:
  #   λ = Σᵢ₌₁ᵏ⁻¹ (αᵢ - αᵢ₊₁) · ωᵢ  +  α_k · ω_k
  # (nodes 1,...,k of A_{n-1})
  R = rank(DT)  # = n-1
  coeffs = zeros(Int, R)
  for i in 1:(k - 1)
    diff = α[i] - α[i + 1]
    if diff != 0
      coeffs[i] = diff  # node i of A_{n-1}
    end
  end
  if α[k] != 0
    coeffs[k] = α[k]    # marked node k (central character = α_k · ω_k)
  end

  CompletelyReducibleBundle(X, coeffs)
end

"""
    schur_functor(X::PartialFlagVariety, partition::NTuple) -> CompletelyReducibleBundle

Schur functor accepting a tuple partition. See also the vector version.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> rank(schur_functor(X, (2, 0)))  # Sym²(U^∨)
3
```
"""
function schur_functor(X::PartialFlagVariety, partition::NTuple{K,<:Integer}) where {K}
  schur_functor(X, collect(Int, partition))
end

# Zero-length partition = structure sheaf
function schur_functor(X::PartialFlagVariety, partition::Tuple{})
  structure_sheaf(X)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Kapranov collection on Grassmannians
# ═══════════════════════════════════════════════════════════════════════════════

"""
    kapranov_bundles_grassmannian(X::PartialFlagVariety) -> Vector{CompletelyReducibleBundle}

All Schur functor bundles ``\\Sigma^\\alpha \\mathcal{U}^\\vee``
on ``\\mathrm{Gr}(k, n)``, ordered by ``|\\alpha| = 0, 1, 2, \\ldots``
(lexicographically within each degree).

These form **Kapranov's full exceptional collection** on ``\\mathrm{Gr}(k, n)``:
the ``\\binom{n}{k}`` bundles indexed by all Young diagrams fitting in a
``k \\times (n-k)`` box.

`X` must be a Type-A generalized Grassmannian.

# References
- M. Kapranov, *On the derived categories of coherent sheaves on some
  homogeneous spaces*, 1988.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> Es = kapranov_bundles_grassmannian(X);

julia> length(Es)
6

julia> euler_characteristic(X)
6
```
"""
function kapranov_bundles_grassmannian(
  X::PartialFlagVariety
)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  DT <: TypeA || throw(ArgumentError(
    "kapranov_bundles_grassmannian requires a Type-A Grassmannian"
  ))
  length(marked) == 1 || throw(
    ArgumentError(
      "kapranov_bundles_grassmannian requires a generalized Grassmannian"
    ),
  )

  k = Int(marked[1])
  n = rank(DT) + 1
  nk = n - k

  # Enumerate all Young diagrams fitting in k × (n-k) box
  partitions = _partitions_in_box(k, nk)

  [schur_functor(X, α) for α in partitions]
end

"""
    _partitions_in_box(k::Int, m::Int) -> Vector{Vector{Int}}

All weakly-decreasing integer sequences ``\\alpha_1 \\geq \\ldots \\geq \\alpha_k \\geq 0``
with ``\\alpha_1 \\leq m``, enumerated in lexicographic order of the tuple.
These are all Young diagrams fitting in a ``k \\times m`` box.
"""
function _partitions_in_box(k::Int, m::Int)
  result = Vector{Int}[]
  _partitions_in_box_helper!(result, Int[], k, m)
  result
end

function _partitions_in_box_helper!(
  result::Vector{Vector{Int}},
  current::Vector{Int},
  remaining::Int,
  max_val::Int,
)
  if remaining == 0
    push!(result, copy(current))
    return nothing
  end
  for v in 0:max_val
    push!(current, v)
    _partitions_in_box_helper!(result, current, remaining - 1, v)
    pop!(current)
  end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Exceptionality predicates on zero loci
#
#  These methods take bundles on the ambient variety X and check
#  exceptionality of their restrictions to the zero locus Z ⊂ X,
#  using the Koszul resolution to compute Ext groups on Z.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    is_exceptional(E::CompletelyReducibleBundle, Z::ZeroLocus) -> Bool

Check whether the restriction ``\\mathcal{E}|_Z`` is an exceptional object on the
zero locus ``Z``.

Computes ``\\operatorname{Ext}^\\vee(\\mathcal{E}|_Z, \\mathcal{E}|_Z) = \\mathrm{H}^\\bullet(Z, \\mathcal{E}^\\vee \\otimes \\mathcal{E}|_Z)``
via the Koszul resolution and checks that ``\\operatorname{Hom} = \\mathbb{C}``
and all higher Ext groups vanish.
"""
function is_exceptional(E::CompletelyReducibleBundle, Z::ZeroLocus)
  endomorphism_bundle = dual(E) ⊗ E
  (H, determined) = cohomology_on_restriction(Z, endomorphism_bundle)
  determined || @warn "cohomology_on_restriction underdetermined for self-Ext"
  H[0] == 1 && all(H[i] == 0 for i in 1:dimension(Z))
end

"""
    is_exceptional_pair(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle,
                        Z::ZeroLocus) -> Bool

Check whether ``(\\mathcal{E}|_Z, \\mathcal{F}|_Z)`` is an exceptional pair on the zero locus ``Z``:
``\\operatorname{RHom}(\\mathcal{F}|_Z, \\mathcal{E}|_Z) = 0``.
"""
function is_exceptional_pair(
  E::CompletelyReducibleBundle,
  F::CompletelyReducibleBundle,
  Z::ZeroLocus,
)
  # Ext^i(F|_Z, E|_Z) = H^i(Z, (F^∨ ⊗ E)|_Z)
  hom_bundle = dual(F) ⊗ E
  (H, determined) = cohomology_on_restriction(Z, hom_bundle)
  determined || @warn "cohomology_on_restriction underdetermined for exceptional pair"
  all(H[i] == 0 for i in 0:dimension(Z))
end

"""
    is_strong_exceptional_pair(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle,
                                Z::ZeroLocus) -> Bool

Check whether ``(\\mathcal{E}|_Z, \\mathcal{F}|_Z)`` is a strong exceptional pair on ``Z``:
``\\operatorname{RHom}(\\mathcal{F}|_Z, \\mathcal{E}|_Z) = 0`` and
``\\operatorname{Ext}^{>0}(\\mathcal{E}|_Z, \\mathcal{F}|_Z) = 0``.
"""
function is_strong_exceptional_pair(
  E::CompletelyReducibleBundle,
  F::CompletelyReducibleBundle,
  Z::ZeroLocus,
)
  is_exceptional_pair(E, F, Z) || return false
  _has_no_positive_exts(E, F, Z)
end

function _has_no_positive_exts(
  E::CompletelyReducibleBundle,
  F::CompletelyReducibleBundle,
  Z::ZeroLocus,
)
  hom_bundle = dual(E) ⊗ F
  (H, determined) = cohomology_on_restriction(Z, hom_bundle)
  determined || @warn "cohomology_on_restriction underdetermined for strong pair"
  all(H[i] == 0 for i in 1:dimension(Z))
end

# Content-based key: the bundles determine the answer, and the defining
# bundle (whose variety is part of its hash) determines the zero locus.
_sequence_status_key(Es::Vector{<:CompletelyReducibleBundle}, Z::ZeroLocus) = (
  collect(CompletelyReducibleBundle, Es), defining_bundle(Z)
)

function _exceptional_sequence_status(
  Es::Vector{<:CompletelyReducibleBundle},
  Z::ZeroLocus,
)
  get!(_SEQUENCE_STATUS_CACHE, _sequence_status_key(Es, Z)) do
    _exceptional_sequence_status_uncached(Es, Z)
  end
end

function _exceptional_sequence_status_uncached(
  Es::Vector{<:CompletelyReducibleBundle},
  Z::ZeroLocus,
)
  n = length(Es)
  all(is_exceptional(Es[i], Z) for i in 1:n) || return (false, false)
  all(is_exceptional_pair(Es[i], Es[j], Z) for i in 1:n for j in (i + 1):n) ||
    return (false, false)
  strong = all(_has_no_positive_exts(Es[i], Es[j], Z) for i in 1:n for j in (i + 1):n)
  (true, strong)
end

"""
    is_exceptional_sequence(Es::Vector{<:CompletelyReducibleBundle},
                            Z::ZeroLocus) -> Bool

Check whether the restrictions ``\\mathcal{E}_i|_Z`` form an exceptional sequence on ``Z``.
"""
function is_exceptional_sequence(
  Es::Vector{<:CompletelyReducibleBundle},
  Z::ZeroLocus,
)
  first(_exceptional_sequence_status(Es, Z))
end

"""
    is_strong_exceptional_sequence(Es::Vector{<:CompletelyReducibleBundle},
                                    Z::ZeroLocus) -> Bool

Check whether the restrictions ``\\mathcal{E}_i|_Z`` form a strong exceptional sequence on ``Z``.
"""
function is_strong_exceptional_sequence(
  Es::Vector{<:CompletelyReducibleBundle},
  Z::ZeroLocus,
)
  last(_exceptional_sequence_status(Es, Z))
end
