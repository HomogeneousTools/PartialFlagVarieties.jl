# ═══════════════════════════════════════════════════════════════════════════════
#  CICY3-1606.04076.jl — Calabi–Yau threefolds in flag varieties
#
#  Classifies Calabi–Yau threefolds as zero loci of sections of
#  completely reducible equivariant bundles on partial flag varieties G/P.
#
#  For Picard rank 1 varieties, we enumerate direct sums of line bundles
#  O(d₁) ⊕ ⋯ ⊕ O(dₖ) with d₁ + ⋯ + dₖ = index(G/P) and k = dim(G/P) - 3.
#
#  Reference: Manivel, "Calabi–Yau threefolds in exceptional flag varieties"
#             (arXiv:1606.04076)
#
#  Usage: julia --project=. examples/CICY3-1606.04076.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Lie

# ═══════════════════════════════════════════════════════════════════════════════
#  Fano index
# ═══════════════════════════════════════════════════════════════════════════════

"""Compute the Fano index of a Picard rank 1 variety G/P."""
function fano_index(X::PartialFlagVariety)
  MDT = marked_type(X)
  Marked = marked_nodes(MDT)
  @assert length(Marked) == 1 "Only Picard rank 1 varieties"
  anticK = PartialFlagVarieties._anticanonical_weight_direct(MDT)
  Int(anticK[Marked[1]])
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Partition enumeration
# ═══════════════════════════════════════════════════════════════════════════════

"""Enumerate all partitions of `n` into exactly `k` parts, each ≥ `min_val`."""
function partitions_into_k(n::Int, k::Int, min_val::Int=1)
  results = Vector{Int}[]
  if k == 0
    n == 0 && push!(results, Int[])
    return results
  end
  if k == 1
    n >= min_val && push!(results, [n])
    return results
  end
  for d in min_val:div(n, k)
    for rest in partitions_into_k(n - d, k - 1, d)
      push!(results, vcat(d, rest))
    end
  end
  results
end

# ═══════════════════════════════════════════════════════════════════════════════
#  CY3 search: complete intersections on Picard-rank-1 varieties
# ═══════════════════════════════════════════════════════════════════════════════

struct CY3Result
  variety_name::String
  dim_ambient::Int
  index::Int
  degrees::Vector{Int}
  chi_top::BigInt
  h11::BigInt
  h21::BigInt
  determined::Bool
end

"""Build the bundle O(d₁) ⊕ ⋯ ⊕ O(dₖ) on X."""
function ci_bundle(X, degrees::Vector{Int})
  E = line_bundle(X, degrees[1])
  for i in 2:length(degrees)
    E = direct_sum(E, line_bundle(X, degrees[i]))
  end
  E
end

"""
Search for CY3 complete intersections on a Picard rank 1 variety.
Returns a vector of CY3Result.
"""
function search_cy3(X, name::String; compute_hodge::Bool=true)
  d = dimension(X)
  d < 4 && return CY3Result[]  # Need dim ≥ 4 for a 3-fold
  idx = fano_index(X)
  codim = d - 3
  parts = partitions_into_k(idx, codim)

  results = CY3Result[]
  for degrees in parts
    E = ci_bundle(X, degrees)

    # Check CY condition
    is_calabi_yau_candidate(E) || continue

    Z = zero_locus(E)

    # Euler characteristic is always exact
    χ_O = euler_characteristic(Z)

    # Compute Hodge numbers if requested and feasible
    h11 = BigInt(-1)
    h21 = BigInt(-1)
    det = false

    if compute_hodge
      try
        h = hodge_numbers(Z)
        h11 = h[2, 2]
        h21 = h[3, 2]
        det = true
      catch e
        # Fall back: compute just χ_top from Ω^p restrictions
      end
    end

    # Topological Euler characteristic: χ_top = 2(h^{1,1} - h^{2,1})
    # (for CY3 with h^{0,0}=h^{3,3}=1, h^{3,0}=h^{0,3}=1, rest 0)
    if h11 >= 0 && h21 >= 0
      chi_top = 2 * (h11 - h21)
    else
      # Compute χ_top = Σ_p (-1)^p χ(Z, Ω^p_X|_Z)
      chi_top = BigInt(0)
      for p in 0:3
        Ωp = exterior_power(cotangent_bundle(X), p)
        chi_top += (-1)^p * euler_characteristic(Z, Ωp)
      end
    end

    push!(results, CY3Result(name, d, idx, degrees, chi_top, h11, h21, det))
  end

  results
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main computation
# ═══════════════════════════════════════════════════════════════════════════════

function main()
  println("=" ^ 72)
  println("  Calabi–Yau threefolds as zero loci in homogeneous varieties G/P")
  println("  Reference: arXiv:1606.04076")
  println("=" ^ 72)
  println()

  all_results = CY3Result[]

  # ── Classical: Projective spaces ──────────────────────────────────────
  println("Searching projective spaces...")
  for n in 4:7
    X = projective_space(n)
    append!(all_results, search_cy3(X, "ℙ$n"))
  end

  # ── Classical: Grassmannians ──────────────────────────────────────────
  println("Searching Grassmannians...")
  for n in 5:7
    for k in 2:div(n, 2)
      d = k * (n - k)
      d >= 4 || continue
      X = Gr(k, n)
      append!(all_results, search_cy3(X, "Gr($k,$n)"))
    end
  end

  # ── Quadrics ──────────────────────────────────────────────────────────
  println("Searching quadrics...")
  for n in 4:6
    X = quadric(n)
    append!(all_results, search_cy3(X, "Q$n"))
  end

  # ── Symplectic/Lagrangian Grassmannians (small) ───────────────────────
  println("Searching symplectic/Lagrangian Grassmannians...")
  for n in [6]
    for k in 2:div(n, 2)
      X = SGr(k, n)
      d = dimension(X)
      d >= 4 || continue
      append!(all_results, search_cy3(X, "SGr($k,$n)"))
    end
  end
  # LGr(3) = SGr(3,6)
  let X = LGr(3)
    append!(all_results, search_cy3(X, "LGr(3)"))
  end

  # ── Orthogonal Grassmannians (small) ──────────────────────────────────
  println("Searching orthogonal Grassmannians...")
  for (k, n) in [(2, 7), (5, 10)]
    X = OGr(k, n)
    d = dimension(X)
    d >= 4 || continue
    append!(all_results, search_cy3(X, "OGr($k,$n)"))
  end

  # ── G₂ varieties ─────────────────────────────────────────────────────
  println("Searching G₂ varieties...")
  let X = partial_flag_variety(TypeG2, 1, "G₂/P₁")
    append!(all_results, search_cy3(X, "G₂/P₁"))
  end
  let X = partial_flag_variety(TypeG2, 2, "G₂/P₂")
    append!(all_results, search_cy3(X, "G₂/P₂"))
  end

  # ── Exceptional: Cayley plane (E₆/P₁) ────────────────────────────────
  println("Searching Cayley plane (dim 16, χ computation only)...")
  let X = cayley_plane()
    idx = fano_index(X)
    d = dimension(X)
    codim = d - 3
    parts = partitions_into_k(idx, codim)
    println("  E₆/P₁: dim=$d, index=$idx, codim=$codim")
    println("  Number of partitions: $(length(parts))")
    # Only compute Euler characteristic (Hodge too expensive for codim 13)
    for degrees in parts
      E = ci_bundle(X, degrees)
      is_calabi_yau_candidate(E) || continue

      Z = zero_locus(E)
      # χ_top via alternating sum of χ(Ω^p|_Z)
      chi_top = BigInt(0)
      for p in 0:3
        Ωp = exterior_power(cotangent_bundle(X), p)
        chi_top += (-1)^p * euler_characteristic(Z, Ωp)
      end
      push!(all_results, CY3Result("OP²", d, idx, degrees, chi_top,
        BigInt(-1), BigInt(-1), false))
    end
  end

  # ── Display results ───────────────────────────────────────────────────
  println()
  println("=" ^ 72)
  println("  Results: $(length(all_results)) CY3 candidates found")
  println("=" ^ 72)
  println()

  if isempty(all_results)
    println("No CY3 candidates found.")
    return
  end

  # Build table data
  names = String[]
  dims = String[]
  indices = String[]
  degree_strs = String[]
  chi_strs = String[]
  h11_strs = String[]
  h21_strs = String[]

  for r in all_results
    push!(names, r.variety_name)
    push!(dims, string(r.dim_ambient))
    push!(indices, string(r.index))
    push!(degree_strs, "(" * join(r.degrees, ",") * ")")
    push!(chi_strs, string(r.chi_top))
    push!(h11_strs, r.h11 >= 0 ? string(r.h11) : "—")
    push!(h21_strs, r.h21 >= 0 ? string(r.h21) : "—")
  end

  data = hcat(names, dims, indices, degree_strs, chi_strs, h11_strs, h21_strs)
  pretty_table(data;
    column_labels=["G/P", "dim", "index", "degrees", "χ_top", "h¹¹", "h²¹"],
    alignment=[:l, :r, :r, :c, :r, :r, :r],
    fit_table_in_display_vertically=false,
  )
end

main()
