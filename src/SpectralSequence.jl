# ═══════════════════════════════════════════════════════════════════════════════
#  SpectralSequence.jl — Spectral Sequences
# ═══════════════════════════════════════════════════════════════════════════════

export SpectralSequence, spectral_sequence, E1_page, isotypical_components, does_E1_degenerate

# ═══════════════════════════════════════════════════════════════════════════════
# Types
# ═══════════════════════════════════════════════════════════════════════════════

"""
    SpectralSequence{T}

A spectral sequence represented by its first page (E₁-page).

# Type Parameters
- `T`: Either `Int` (for integer multiplicities) or `WeylCharacter` (for character-valued entries)

# Fields
- `E1::Dict{Tuple{Int, Int}, T}`: The E₁-page as a dictionary mapping bidegree positions `(p, q)` 
  to entries of type `T`. Position `(p, q)` represents total degree `p + q` with vertical degree `q`.
"""

struct SpectralSequence{T <: Union{Int, WeylCharacter}}
  E1::Dict{Tuple{Int, Int}, T}
  function SpectralSequence{T}(E1::Dict{Tuple{Int, Int}, T}) where {T <: Union{Int, WeylCharacter}}
    new{T}(E1)
  end
end

"""
    E1_page(S::SpectralSequence{T}) -> Dict{Tuple{Int, Int}, T}

Extract the E₁-page from a spectral sequence.
"""

E1_page(S::SpectralSequence{T}) where {T} = S.E1
# ═══════════════════════════════════════════════════════════════════════════════
#  Constructors
# ═══════════════════════════════════════════════════════════════════════════════

"""
    spectral_sequence(F::FilteredBundle) -> SpectralSequence{WeylCharacter}

Given a filtered bundle with graded pieces, this function computes the E₁-page by calculating
the sheaf cohomology of each graded piece using the Borel–Weil–Bott theorem.

# Details
The E₁-page entry at position `(p, q)` contains the sheaf cohomology in degree `p+q` of the 
`q`-th graded piece of the filtration. In other words, the spectral sequence uses the convention where
bidegree is `(-q + i, q)` for filtration index `q` and cohomology degree `i`.
"""
function spectral_sequence(F::FilteredBundle)
  n = n_filtration_steps(F)
  bundles = reverse(graded_pieces(F)) # Reverse to match spectral sequence convention
  d = dimension(variety(F))
  DT = dynkin_type(variety(F))
  R = rank(variety(F))
  zero_character = WeylCharacter(DT)
  E = Dict{Tuple{Int, Int}, WeylCharacter{DT,R}}()
  for q in 0:n-1
    H = cohomology(bundles[q+1])
    for i in 0:d
      H[i] != zero_character && (E[(-q+i,q)] = H[i])
    end
  end
  return SpectralSequence{WeylCharacter{DT,R}}(E)
end

"""
    isotypical_components(S::SpectralSequence{WeylCharacter}) -> Dict{WeightLatticeElem, SpectralSequence{Int}}

Decompose a spectral sequence into its isotypical components.

For each weight in the weight lattice that appears in some E₁-page entry, this function
extracts the multiplicity of that weight across the entire spectral sequence, yielding a
new spectral sequence with integer entries. This decomposition respects the bidegree structure.

# Details
Each isotypical component is a spectral sequence where the entry at position `(p, q)` is the
multiplicity of the weight in the character at that position (or 0 if the weight does not appear).
"""
function isotypical_components(S::SpectralSequence{WeylCharacter{DT,R}}) where {DT, R}
  # Accumulate components in a single pass through E1
  iso_dict = Dict{WeightLatticeElem{DT,R}, Dict{Tuple{Int,Int}, Int}}()
  
  for (pos, char) in S.E1
    for (weight, mult) in char.terms
      if !haskey(iso_dict, weight)
        iso_dict[weight] = Dict{Tuple{Int,Int}, Int}()
      end
      iso_dict[weight][pos] = mult
    end
  end
  
  # Convert accumulated dicts to SpectralSequence objects
  return Dict(weight => SpectralSequence{Int}(E_iso) 
              for (weight, E_iso) in iso_dict)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Degeneracy Check
# ═══════════════════════════════════════════════════════════════════════════════

"""
    does_E1_degenerate(S::SpectralSequence{Int}) -> Bool

Test whether the E₁-page of an integer-valued spectral sequence degenerates.

# Details
The function checks for the existence of pairs of positions `(p₁, q₁)` and `(p₂, q₂)` such that:
- Total degrees satisfy `p₁ + q₁ = p₂ + q₂ - 1`
- Vertical degrees satisfy `q₁ < q₂`
"""
function does_E1_degenerate(S::SpectralSequence{Int})
  E = E1_page(S)
  return !any(pos_1[1] + pos_1[2] == pos_2[1] + pos_2[2] - 1 && pos_1[2] < pos_2[2] 
              for pos_1 in keys(E) 
              for pos_2 in keys(E))
end

"""
    does_E1_degenerate(S::SpectralSequence{WeylCharacter}) -> Bool

Test whether the E₁-page of a character-valued spectral sequence degenerates at E₁.

# Details
The function first decomposes the spectral sequence into isotypical components (one for each
weight in the weight lattice), then checks whether each component degenerates using the
integer-valued degeneracy test. If any component has a nonzero differential, returns `false`.

# Examples
```jldoctest

julia> using PartialFlagVarieties

julia> T = tangent_bundle_filtration(SGr(4, 10))

julia> S_1 = spectral_sequence(exterior_power(T,2))

julia> S_2 = spectral_sequence(exterior_power(T,5))

julia> does_E1_degenerate(S_1)
true

julia> does_E1_degenerate(S_2)
false
```
"""
function does_E1_degenerate(S::SpectralSequence{WeylCharacter{DT,R}}) where {DT, R}
  iso = isotypical_components(S)
  return all(does_E1_degenerate, values(iso))
end