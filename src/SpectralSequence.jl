# ═══════════════════════════════════════════════════════════════════════════════
#  SpectralSequence.jl — Spectral Sequences
# ═══════════════════════════════════════════════════════════════════════════════

export SpectralSequence,
  spectral_sequence, E1_page, isotypical_components, does_E1_degenerate

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

struct SpectralSequence{T<:Union{Int,WeylCharacter}}
  E1::Dict{Tuple{Int,Int},T}
  function SpectralSequence{T}(
    E1::Dict{Tuple{Int,Int},T}
  ) where {T<:Union{Int,WeylCharacter}}
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
  E = Dict{Tuple{Int,Int},WeylCharacter{DT,R}}()
  for q in 0:(n - 1)
    H = cohomology(bundles[q + 1])
    for i in 0:d
      H[i] != zero_character && (E[(-q+i, q)] = H[i])
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
function isotypical_components(S::SpectralSequence{WeylCharacter{DT,R}}) where {DT,R}
  # Accumulate components in a single pass through E1
  iso_dict = Dict{WeightLatticeElem{DT,R},Dict{Tuple{Int,Int},Int}}()

  for (pos, char) in S.E1
    for (weight, mult) in char.terms
      if !haskey(iso_dict, weight)
        iso_dict[weight] = Dict{Tuple{Int,Int},Int}()
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
  return !any(
    pos_1[1] + pos_1[2] == pos_2[1] + pos_2[2] - 1 && pos_1[2] < pos_2[2]
    for pos_1 in keys(E)
    for pos_2 in keys(E)
  )
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
function does_E1_degenerate(S::SpectralSequence{WeylCharacter{DT,R}}) where {DT,R}
  iso = isotypical_components(S)
  return all(does_E1_degenerate, values(iso))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Higher Cohomology Check (Experimental)
# ═══════════════════════════════════════════════════════════════════════════════
"""
    has_higher_cohomology(S::SpectralSequence{Int}) -> Bool

Test whether a filtered bundle has nonzero higher cohomology.

# Details:

We inspect the E₁-page to detect whether higher cohomology (at positions where `p + q ≠ 0`) is nonzero.

# TODO

Examples, better name.

# WARNING

If the function returns `false`, it does not guarantee that the filtered bundle has no higher cohomology.
"""
function has_higher_cohomology(S::SpectralSequence{Int})
  return any(
    does_E1_degenerate(S) &&
    any(pos[1] + pos[2] != 0 for pos in keys(E1_page(S))),
  )
end

"""
    has_higher_cohomology_experimental(S::SpectralSequence{Int}) -> Bool

Test whether a filtered bundle has nonzero higher cohomology by analyzing
the E₁-page of the spectral sequence without requiring full degeneracy.

# Details:

If `E1[p,q]` is nonzero for some `p+q≠0`, checks all possible maps to `E1[p,q]` that could 
potentially kill the cohomology. If the total dimension of these potentially killing maps is less 
than the dimension of `E1[p,q]`, then we can conclude that some cohomology survives to E∞, and 
return `true`. If we find no such evidence, we return `false`.

# TODO

Better documentation, change name, examples?
"""
function has_higher_cohomology_experimental(S::SpectralSequence{Int})
  for pos1 in keys(E1_page(S))
    total_degree = pos1[1] + pos1[2]

    # Skip H⁰ 
    if total_degree == 0
      continue
    end

    # Sum dimensions of positions that could potentially kill cohomology at pos1
    killing_dimension = sum(
      (
        S.E1[pos2] for pos2 in keys(E1_page(S))
        if (pos2[1] + pos2[2] > total_degree && pos1[2] < pos2[2]) ||
        (pos2[1] + pos2[2] < total_degree && pos1[2] > pos2[2])
      );
      init=0,
    )

    # If killing dimension < dimension at pos1, cohomology survives to E∞
    if killing_dimension < S.E1[pos1]
      return true
    end
  end
  return false
end
"""
    has_higher_cohomology_experimental_2(S::SpectralSequence{Int}) -> Bool

Test whether a filtered bundle has nonzero higher cohomology by analyzing
the E₁-page of the spectral sequence without requiring full degeneracy.

We want to check if it is possible to kill all cohomology above the diagonal, which are exactly the 
positions associated with higher cohomology. If we cannot kill all cohomology above the diagonal,
then we can conclude that some cohomology survives to E∞, and return `true`.

# Details:

1. If the page is concentrated on the diagonal, return `false`.
2. If there is any entry (not on the diagonal) that has no possible differential maps in or out, return `true`.
3. Look for positions where there is only one differential map that might kill that cohomology.
   We create a new spectral sequence where we killed the cohomology at that position, go back to step 1,
   and repeat.
4. If there are no positions with only zero or one possible killing map, return `false`.

# TODO

Better documentation, change name, examples?
"""
function has_higher_cohomology_experimental_2(S::SpectralSequence{Int})
  E1 = E1_page(S)
  #If there only diagonal entries, return false
  if all(x[1]+x[2] == 0 for x in keys(E1))
    return false
  end
  for x in keys(E1)
    total_degree = x[1] + x[2]
    # We skip diagonal entries.
    if (total_degree==0)
      continue
    end
    # We only consider pairs of entries that could potentially kill each other.
    killers = Vector{Tuple{Int,Int}}()
    for y in keys(E1)
      if (
        (total_degree == y[1]+y[2]-1 && x[2] < y[2]) ||
        (total_degree == y[1]+y[2]+1 && x[2] > y[2])
      )
        push!(killers, y)
      end
    end
    if (length(killers) == 0)
      # If there are no killers, then the cohomology class at x survives to E∞, so we return true.
      return true
    end
    if (length(killers) == 1)
      killer = killers[1]
      new_E1 = copy(E1)
      if E1[killer] > E1[x]
        new_E1[killer] -= E1[x]
        delete!(new_E1, x)
      elseif E1[killer] < E1[x]
        new_E1[x] -= E1[killer]
        delete!(new_E1, killer)
      else
        delete!(new_E1, x)
        delete!(new_E1, killer)
      end
      return has_higher_cohomology_experimental_2(SpectralSequence{Int}(new_E1))
    end
  end
  #If we ever leave the loop, it means that there is never a unique killer, so
  #we cannot conclude that any cohomology class survives, so we return false.
  return false
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Cohomology of FilteredBundle via spectral sequence
# TODO: most of it IA generated, check carefully
# ═══════════════════════════════════════════════════════════════════════════════

"""
Collect the diagonal sums for an integer-valued spectral sequence:
`D[i] = Σ_{p+q=i} E1[(p,q)]`.
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
Return `true` if there exist positions `pos_1` on total-degree diagonal `i`
and `pos_2` on diagonal `i+1` in `S` satisfying `pos_1[2] < pos_2[2]`.
This is the per-diagonal analogue of the `does_E1_degenerate` check.
"""
function _has_potential_diff(S::SpectralSequence{Int}, i::Int)
  E = E1_page(S)
  q_i = [pos[2] for pos in keys(E) if pos[1] + pos[2] == i]
  q_j = [pos[2] for pos in keys(E) if pos[1] + pos[2] == i + 1]
  isempty(q_i) && return false
  isempty(q_j) && return false
  return minimum(q_i) < maximum(q_j)
end

"""
Group symbolic cohomology entries by their coefficient signature and re-express
each group using a single fresh variable.

Two entries `h^i = c_i + P` and `h^j = c_j + P` that share the same symbolic
part `P` (identical `coeffs` dict) are replaced by `y + (c_i - c_min)` and
`y + (c_j - c_min)`, where `y` is a new variable and `c_min` is the smallest
constant in the group.  Singleton groups are similarly flattened to `y + 0`.

After this call every non-determined entry contains exactly one variable with
coefficient +1.  A subsequent `_renumber_variables!` compacts the IDs to 0,1,…
"""
function _reduce_cohomology_entries!(entries::Vector{AffineExpr}, var_counter::Ref{Int})
  sig_to_indices = Dict{Vector{Pair{Int,BigInt}},Vector{Int}}()
  for (idx, e) in enumerate(entries)
    isempty(e.coeffs) && continue
    sig = sort!(collect(e.coeffs), by=first)
    push!(get!(sig_to_indices, sig, Int[]), idx)
  end
  for (_, group) in sig_to_indices
    pivot_const = minimum(entries[i].constant for i in group)
    new_id = var_counter[]
    var_counter[] += 1
    for idx in group
      delta = entries[idx].constant - pivot_const
      entries[idx] = AffineExpr(delta, Dict{Int,BigInt}(new_id => BigInt(1)))
    end
  end
end

"""
    cohomology(F::FilteredBundle) -> (Cohomology{AffineExpr}, Bool)

Compute the cohomology of a filtered bundle via its spectral sequence.

For each isotypical component of the E₁-page:
- If consecutive non-zero diagonals `i` and `i+1` admit a potential differential
  (i.e., some position on diagonal `i` has strictly smaller `q`-coordinate than
  some position on diagonal `i+1`), a symbolic variable is introduced and
  subtracted from both the `h^i` and `h^{i+1}` contributions of that component.
- The total `h^i` is the sum over all isotypical components of
  `(diagonal_sum_i - x_{i-1,i} - x_{i,i+1}) * dim(V_w)`.

After assembling the raw expressions, entries that share the same symbolic
structure are reduced to a single free parameter (see `_reduce_cohomology_entries!`).
This means if `h^i` and `h^j` differ only by a constant, the output uses one
variable: `h^i = x_0`, `h^j = x_0 + (c_j - c_i)`.

Returns `(H, true)` when every isotypical component degenerates (no variables
are needed, so `H` is fully determined), and `(H, false)` otherwise.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> T = tangent_bundle_filtration(SGr(4, 10))

julia> H, det = cohomology(exterior_power(T, 2));

julia> det
true
```
"""
function cohomology(F::FilteredBundle)
  S = spectral_sequence(F)
  d = dimension(variety(F))
  iso = isotypical_components(S)
  var_counter = Ref(0)
  entries = AffineExpr[AffineExpr(0) for _ in 0:d]
  all_determined = true

  for (weight, S_iso) in iso
    dim_w = BigInt(degree(weight))
    D = _diagonal_sums(S_iso)
    total_degs = sort(collect(keys(D)))

    # Assign a symbolic variable to each consecutive pair of diagonals
    # that admits a potential differential.
    pair_vars = Dict{Tuple{Int,Int},AffineExpr}()
    for k in 1:(length(total_degs) - 1)
      i, j = total_degs[k], total_degs[k + 1]
      j == i + 1 || continue
      _has_potential_diff(S_iso, i) || continue
      pair_vars[(i, j)] = symbolic_variable(var_counter[])
      var_counter[] += 1
      all_determined = false
    end

    # Accumulate contributions into the cohomology entries.
    for (i, d_i) in D
      0 <= i <= d || continue
      contribution = AffineExpr(d_i)
      haskey(pair_vars, (i - 1, i)) && (contribution -= pair_vars[(i - 1, i)])
      haskey(pair_vars, (i, i + 1)) && (contribution -= pair_vars[(i, i + 1)])
      entries[i + 1] += dim_w * contribution
    end
  end

  # Reduce: entries sharing the same symbolic structure collapse to one variable.
  if !all_determined
    _reduce_cohomology_entries!(entries, var_counter)
    mat = reshape(entries, 1, length(entries))
    _renumber_variables!(mat)
  end

  (Cohomology{AffineExpr}(entries, d), all_determined)
end