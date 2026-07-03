# ═══════════════════════════════════════════════════════════════════════════════
#  SpectralSequence.jl — the spectral sequence of a filtered bundle
#
#  A FilteredBundle F on X = G/P determines a spectral sequence
#
#      E₁ = H^*(X, gr F)  ⟹  H^*(X, F)
#
#  computed from the graded pieces via Borel–Weil–Bott.  We store the E₁
#  page at positions (p, q) where q counts the graded pieces from the top
#  of the filtration and p + q is the cohomological degree; every
#  differential then goes from (p, q) to (p + 1 - r, q + r) with r ≥ 1,
#  i.e. raises the total degree by 1 and strictly raises q.
#
#  All differentials are G-equivariant, so the spectral sequence splits
#  into isotypical components indexed by the highest weights appearing on
#  the E₁ page.  Degeneracy tests and dimension bookkeeping are done per
#  component, which is much sharper than counting dimensions alone.
# ═══════════════════════════════════════════════════════════════════════════════

export SpectralSequence,
  spectral_sequence, E1_page, isotypical_components, does_E1_degenerate

# ═══════════════════════════════════════════════════════════════════════════════
#  Type
# ═══════════════════════════════════════════════════════════════════════════════

"""
    SpectralSequence{T}

A spectral sequence represented by its E₁ page, stored as a dictionary
mapping positions `(p, q)` to entries of type `T`.  The entry at `(p, q)`
sits in total degree `p + q`; differentials raise the total degree by 1 and
strictly raise `q`.

For the spectral sequence of a filtered bundle the support is **not** the
first quadrant but the horizontal band

```math
0 \\le q \\le s - 1, \\qquad 0 \\le p + q \\le \\dim X,
```

where ``s`` is the number of graded pieces: ``q`` indexes the pieces from
the top of the filtration, the total degree is a cohomological degree, and
``p = (p+q) - q`` may therefore be negative.

The type parameter is either `WeylCharacter` (character-valued entries) or
`Int` (multiplicities, as in an isotypical component).
"""
struct SpectralSequence{T<:Union{Int,WeylCharacter}}
  E1::Dict{Tuple{Int,Int},T}
end

"""
    E1_page(S::SpectralSequence{T}) -> Dict{Tuple{Int, Int}, T}

The E₁ page of the spectral sequence.
"""
E1_page(S::SpectralSequence) = S.E1

# ═══════════════════════════════════════════════════════════════════════════════
#  Constructor from a filtered bundle
# ═══════════════════════════════════════════════════════════════════════════════

"""
    spectral_sequence(F::FilteredBundle) -> SpectralSequence{WeylCharacter}

The spectral sequence ``E_1 = \\mathrm{H}^\\bullet(X, \\mathrm{gr}\\, F) \\Rightarrow \\mathrm{H}^\\bullet(X, F)``
of a filtered bundle, with the E₁ page computed by Borel–Weil–Bott.

The cohomology of the `q`-th graded piece **counted from the top** of the
filtration is placed in positions of second coordinate `q` (so `q = 0` is
the top quotient), with the cohomological degree as total degree.  With
this convention every differential raises the total degree by 1 and
strictly raises `q`: for a two-step filtration ``0 \\to F_1 \\to F \\to
\\mathrm{gr}_2 \\to 0`` the connecting map ``\\mathrm{H}^i(\\mathrm{gr}_2) \\to
\\mathrm{H}^{i+1}(F_1)`` goes from `q = 0` to `q = 1`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> F = filtered_tangent_bundle(SGr(2, 6));

julia> S = spectral_sequence(F);

julia> E1_page(S)  # H⁰(T) is the adjoint representation, from the top piece
Dict{Tuple{Int64, Int64}, WeylCharacter{TypeC{3}, 3}} with 1 entry:
  (0, 0) => C3(2, 0, 0)
```
"""
function spectral_sequence(F::FilteredBundle)
  DT = dynkin_type(variety(F))
  R = rank(variety(F))
  d = dimension(variety(F))
  E = Dict{Tuple{Int,Int},WeylCharacter{DT,R}}()
  for (k, piece) in enumerate(reverse(graded_pieces(F)))
    q = k - 1
    H = cohomology(piece)
    for i in 0:d
      isempty(H[i].terms) || (E[(i - q, q)] = H[i])
    end
  end
  SpectralSequence(E)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Isotypical decomposition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    isotypical_components(S::SpectralSequence{WeylCharacter})
      -> Dict{WeightLatticeElem, SpectralSequence{Int}}

Decompose a character-valued spectral sequence into its isotypical
components: for each highest weight appearing on the E₁ page, the
integer-valued spectral sequence of its multiplicities.

The differentials are ``G``-equivariant, so they preserve this
decomposition; any statement proved componentwise (degeneracy, surviving
dimensions) holds for the full spectral sequence.
"""
function isotypical_components(S::SpectralSequence{WeylCharacter{DT,R}}) where {DT,R}
  components = Dict{WeightLatticeElem{DT,R},Dict{Tuple{Int,Int},Int}}()
  for (pos, char) in S.E1
    for (weight, mult) in char.terms
      get!(components, weight, Dict{Tuple{Int,Int},Int}())[pos] = mult
    end
  end
  Dict(weight => SpectralSequence(E) for (weight, E) in components)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Degeneracy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    does_E1_degenerate(S::SpectralSequence) -> Bool

Test whether the spectral sequence has no room for a nonzero differential,
i.e. no pair of E₁ entries in adjacent total degrees with strictly
increasing `q`.

This is a sufficient condition for degeneracy at E₁: when it holds,
``E_1 = E_\\infty``.  When it fails the spectral sequence may still
degenerate for other reasons, so `false` means "unknown".

For a character-valued spectral sequence the test is applied to each
isotypical component separately, which is sharper: a differential must be
equivariant, so entries of different weights cannot interact.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = SGr(2, 6); Omega = dual(filtered_tangent_bundle(X));

julia> does_E1_degenerate(spectral_sequence(exterior_power(Omega, 3)))
true

julia> W = tensor_product(exterior_power(Omega, 3), line_bundle(X, -1));

julia> does_E1_degenerate(spectral_sequence(W))
false
```
"""
function does_E1_degenerate(S::SpectralSequence{Int})
  diagonals = Set(p[1] + p[2] for p in keys(E1_page(S)))
  !any(_has_potential_diff(S, i) for i in diagonals)
end

function does_E1_degenerate(S::SpectralSequence{<:WeylCharacter})
  all(does_E1_degenerate, values(isotypical_components(S)))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Cohomology of a filtered bundle
# ═══════════════════════════════════════════════════════════════════════════════

"""
Total dimension of an integer-valued spectral sequence on each total-degree
diagonal: `D[i] = Σ_{p+q=i} E1[(p,q)]`.
"""
function _diagonal_sums(S::SpectralSequence{Int})
  D = Dict{Int,Int}()
  for (pos, mult) in S.E1
    i = pos[1] + pos[2]
    D[i] = get(D, i, 0) + mult
  end
  D
end

"""
Whether some differential can connect total degree `i` to total degree
`i + 1` in `S`: true when an entry on diagonal `i` has strictly smaller `q`
than an entry on diagonal `i + 1`.
"""
function _has_potential_diff(S::SpectralSequence{Int}, i::Int)
  E = E1_page(S)
  q_from = [pos[2] for pos in keys(E) if pos[1] + pos[2] == i]
  q_to = [pos[2] for pos in keys(E) if pos[1] + pos[2] == i + 1]
  !isempty(q_from) && !isempty(q_to) && minimum(q_from) < maximum(q_to)
end

"""
    _cohomology_filtered(F::FilteredBundle, var_counter) -> Vector{AffineExpr}

``\\mathrm{H}^\\bullet(X, F)`` via the spectral sequence of the filtration, as a vector of
affine expressions indexed by cohomological degree `0:dim(X)`.

For each isotypical component, differentials can only cancel E₁ classes on
adjacent total-degree diagonals, in equal amounts on both sides.  One
symbolic variable is introduced per pair of adjacent diagonals that admits
a potential differential; it stands for the **total dimension** cancelled
between the two diagonals (the weight multiplicity times the dimension of
the irreducible), so it enters the entries with coefficient ±1 and integer
elimination in `_apply_equation!` always applies.  Components whose
spectral sequence visibly degenerates contribute exact dimensions.
"""
function _cohomology_filtered(F::FilteredBundle, var_counter::Ref{Int})
  d = dimension(variety(F))

  # A one-step filtration has nothing to degenerate: E₁ = E_∞ exactly.
  if n_filtration_steps(F) <= 1
    H = dimensions(total_bundle(F))
    return AffineExpr[AffineExpr(H[i]) for i in 0:d]
  end

  entries = AffineExpr[AffineExpr(0) for _ in 0:d]
  for (weight, S) in isotypical_components(spectral_sequence(F))
    dim_w = BigInt(degree(weight))
    D = _diagonal_sums(S)

    # cancelled[i] = total dimension killed between diagonals i and i + 1
    cancelled = Dict{Int,AffineExpr}()
    for i in keys(D)
      if haskey(D, i + 1) && _has_potential_diff(S, i)
        cancelled[i] = _fresh_variable(var_counter)
      end
    end

    for (i, D_i) in D
      0 <= i <= d || continue
      e = AffineExpr(dim_w * D_i)
      haskey(cancelled, i - 1) && (e -= cancelled[i - 1])
      haskey(cancelled, i) && (e -= cancelled[i])
      entries[i + 1] += e
    end
  end
  entries
end

"""
    cohomology(F::FilteredBundle) -> Cohomology{AffineExpr}

Compute ``\\mathrm{H}^\\bullet(X, F)`` via the spectral sequence of the filtration.

Entries are affine expressions: exact integers where the spectral sequence
visibly degenerates (checked per isotypical component), and expressions in
symbolic variables ``x_0, x_1, \\ldots`` — one per pair of adjacent
total-degree diagonals of an isotypical component that admits a potential
differential — otherwise.  Use `is_determined` to check individual entries.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = SGr(2, 6); Omega = dual(filtered_tangent_bundle(X));

julia> cohomology(exterior_power(Omega, 2))
H² = 2

julia> cohomology(tensor_product(exterior_power(Omega, 3), line_bundle(X, -1)))
H⁴ = 2 - x_0
H⁵ = 1 - x_0
```
"""
function cohomology(F::FilteredBundle)
  entries = _cohomology_filtered(F, Ref(0))
  _renumber_variables!(entries)
  Cohomology{AffineExpr}(entries, dimension(variety(F)))
end
