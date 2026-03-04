# ═══════════════════════════════════════════════════════════════════════════════
#  CICY3-1606.04076.jl — Calabi–Yau threefolds as zero loci
#
#  Searches for Calabi–Yau threefolds arising as zero loci of sections
#  of completely reducible equivariant vector bundles on:
#    • Projective spaces ℙ^n  (n = 4,…,7)
#    • Grassmannians Gr(2,n)  (n = 5,…,7)
#    • Quadrics Q_n           (n = 5,…,7)
#    • Lagrangian Grassmannians LGr(n,2n)  (n = 3,4)
#    • Orthogonal Grassmannians OGr(n,2n)  (n = 5)
#    • G₂-Grassmannian G₂/P₁ (5-dimensional quadric)
#    • Cayley plane OP²
#
#  The classification recovers results from:
#    Benedetti–Fatighenti–Mongardi–Ruiperez, "Calabi-Yau threefolds in
#    flag varieties" (arXiv:1606.04076)
#
#  For each ambient variety X, we enumerate completely reducible bundles E
#  with rk(E) = dim(X) - 3 (so Z(s) is a threefold) and check the
#  Calabi–Yau condition det(E) ≅ ω_X^{-1}.
#
#  Usage: julia --project=. examples/CICY3-1606.04076.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Lie

# ═══════════════════════════════════════════════════════════════════════════════
#  Enumeration of completely reducible bundles
# ═══════════════════════════════════════════════════════════════════════════════

"""
    enumerate_cr_bundles(X, irreps, target_rank, target_det)

Enumerate completely reducible bundles on `X` with given `target_rank`
and `target_det` (determinant central character), using summands from
`irreps` (a list of (IrrepLevi, rank, det_central) triples).

Returns a Vector of CompletelyReducibleBundle.
"""
function enumerate_cr_bundles(
  X::PartialFlagVariety{MDT},
  irreps::Vector{<:Tuple},
  target_rank::Int,
  target_det::Vector{Rational{Int}},
) where {MDT}
  results = CompletelyReducibleBundle{MDT}[]

  function backtrack(idx, remaining_rank, remaining_det, summands)
    if remaining_rank == 0 && all(remaining_det .== 0)
      # Wrap IrrepLevi summands as CompletelyReducibleBundle and combine
      bundles = [CompletelyReducibleBundle{MDT}(X, [s]) for s in summands]
      E = reduce(direct_sum, bundles)
      push!(results, E)
      return
    end
    if remaining_rank <= 0 || idx > length(irreps)
      return
    end

    irr, rk, det_c = irreps[idx]
    rk_int = Int(rk)

    # Maximum multiplicity of this summand
    max_mult = remaining_rank ÷ rk_int

    for m in 0:max_mult
      new_rank = remaining_rank - m * rk_int
      new_det = remaining_det .- m .* det_c

      new_summands = copy(summands)
      for _ in 1:m
        push!(new_summands, irr)
      end

      backtrack(idx + 1, new_rank, new_det, new_summands)
    end
  end

  backtrack(1, target_rank, target_det, IrrepLevi{MDT}[])
  results
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Candidate irreducible summands for each variety type
# ═══════════════════════════════════════════════════════════════════════════════

"""
Get irreducible representations to use as summands, together with their
rank and determinant central character.
"""
function get_candidate_irreps(X::PartialFlagVariety{MDT}; max_rank::Int=0) where {MDT}
  DT = PartialFlagVarieties._ambient_type(MDT)
  R = Lie.rank(DT)
  Marked = marked_nodes(MDT)

  candidates = Tuple{IrrepLevi{MDT},BigInt,Vector{Rational{Int}}}[]

  # For each marked node m, the line bundle O(1) has weight ω_m
  # We consider line bundles O(d) for d = 1, 2, 3, ...
  # and low-dimensional irreps

  # Strategy: build a list of candidate irreps:
  # 1) Line bundles O(d) for each marked node
  # 2) Small exterior powers of universal quotient/subbundle duals

  # d_max controls the maximum degree of line bundles O(d) we try.
  # Must be large enough to capture all relevant summands.
  # For line bundles, total degree = sum(d_i) must match anticanonical degree.
  d_max = max_rank > 0 ? max(max_rank, 8) : 8

  for m in Marked
    for d in 1:d_max
      omega = zeros(Int, R)
      omega[m] = d
      λ = WeightLatticeElem(DT, omega)
      irr = IrrepLevi(MDT, λ)

      rk = rank_bundle(CompletelyReducibleBundle{MDT}(X, [irr]))
      det_c = collect(PartialFlagVarieties._determinant_central(
        CompletelyReducibleBundle{MDT}(X, [irr])))

      push!(candidates, (irr, rk, det_c))
    end
  end

  # Add small exterior powers for Grassmannians (2 marked nodes = full flag,
  # but for Gr(k,n) we have 1 marked node)
  if length(Marked) == 1
    m = Marked[1]

    # For type A: ω_i gives ∧^i(S*) if i ≤ k, or ∧^i(S*) ⊗ det
    # We just enumerate small fundamental representations
    for i in 1:R
      if i != m  # ω_m is already covered as O(1)
        omega = zeros(Int, R)
        omega[i] = 1
        λ = WeightLatticeElem(DT, omega)
        irr = IrrepLevi(MDT, λ)

        rk = rank_bundle(CompletelyReducibleBundle{MDT}(X, [irr]))
        det_c = collect(PartialFlagVarieties._determinant_central(
          CompletelyReducibleBundle{MDT}(X, [irr])))

        push!(candidates, (irr, rk, det_c))
      end
    end
  end

  # Filter: only keep candidates with rank ≤ target
  if max_rank > 0
    filter!(c -> c[2] <= max_rank, candidates)
  end

  # Remove duplicates (same irrep)
  unique!(c -> c[1], candidates)

  candidates
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Result type
# ═══════════════════════════════════════════════════════════════════════════════

struct CY3Result
  ambient::String
  dim_ambient::Int
  bundle_desc::String
  rank_E::Int
  chi_O::BigInt
  h11::BigInt
  h21::BigInt
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Search on a given variety
# ═══════════════════════════════════════════════════════════════════════════════

function search_cy3(X::PartialFlagVariety{MDT}, name::String, results::Vector{CY3Result};
  max_summand_rank::Int=0) where {MDT}
  d = dimension(X)
  codim = d - 3
  codim >= 1 || return

  print("  Searching $name (dim=$d, codim=$codim)... ")
  flush(stdout)

  target_rank = Int(codim)

  # Compute anticanonical central character = target determinant
  antican = collect(PartialFlagVarieties._anticanonical_central(MDT))

  # Get candidate irreps
  irreps = get_candidate_irreps(X; max_rank=max_summand_rank > 0 ? max_summand_rank : target_rank)

  # Enumerate bundles with correct rank and CY determinant condition
  bundles = enumerate_cr_bundles(X, irreps, target_rank, antican)

  println("found $(length(bundles)) CY candidate(s)")

  for E in bundles
    Z = zero_locus(E)
    chi_O = euler_characteristic(Z)

    # Attempt Hodge numbers
    h11 = BigInt(-1)
    h21 = BigInt(-1)
    try
      h = hodge_numbers(Z)
      h11 = h[2, 2]  # h^{1,1}
      h21 = h[2, 3]  # h^{1,2} = h^{2,1}
    catch e
      println("    (Hodge failed for $E: $(sprint(showerror, e)))")
    end

    # Build bundle description from summands
    desc = join([string(s) for s in E.components], " ⊕ ")

    push!(results, CY3Result(name, d, desc, Int(rank_bundle(E)), chi_O, h11, h21))
  end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main computation
# ═══════════════════════════════════════════════════════════════════════════════

function main()
  println("=" ^ 72)
  println("  Calabi–Yau threefolds as zero loci in flag varieties")
  println("  Following Benedetti–Fatighenti–Mongardi–Ruiperez (1606.04076)")
  println("=" ^ 72)
  println()

  results = CY3Result[]

  # ── Projective spaces ──────────────────────────────────────────────
  println("Projective spaces:")
  for n in 4:7
    X = projective_space(n)
    search_cy3(X, "ℙ$n", results)
  end
  println()

  # ── Grassmannians ──────────────────────────────────────────────────
  println("Grassmannians:")
  for (k, n) in [(2, 5), (2, 6), (2, 7)]
    X = Gr(k, n)
    search_cy3(X, "Gr($k,$n)", results; max_summand_rank=Int(dimension(X) - 3))
  end
  println()

  # ── Quadrics ───────────────────────────────────────────────────────
  println("Quadrics:")
  for n in 5:7
    X = quadric(n)
    search_cy3(X, "Q$n", results)
  end
  println()

  # ── Lagrangian Grassmannians ───────────────────────────────────────
  println("Lagrangian/Symplectic Grassmannians:")
  for n in 3:4
    X = LGr(n)
    search_cy3(X, "LGr($n,$(2n))", results)
  end
  println()

  # ── Orthogonal Grassmannians ───────────────────────────────────────
  println("Orthogonal Grassmannians:")
  X = OGr(5, 10)
  search_cy3(X, "OGr(5,10)", results)
  println()

  # ── Cayley plane ───────────────────────────────────────────────────
  println("Exceptional varieties:")
  try
    X = cayley_plane()
    search_cy3(X, "OP²", results)
  catch e
    println("  OP² skipped: $(sprint(showerror, e))")
  end
  println()

  # ── Display results ────────────────────────────────────────────────
  println("=" ^ 72)
  println("  Results: $(length(results)) Calabi–Yau threefold(s) found")
  println("=" ^ 72)
  println()

  if isempty(results)
    println("  No CY3 candidates found.")
    return
  end

  ambients = [r.ambient for r in results]
  descs = [r.bundle_desc for r in results]
  ranks = [string(r.rank_E) for r in results]
  chis = [string(r.chi_O) for r in results]
  h11s = [r.h11 >= 0 ? string(r.h11) : "—" for r in results]
  h21s = [r.h21 >= 0 ? string(r.h21) : "—" for r in results]

  data = hcat(ambients, descs, ranks, chis, h11s, h21s)
  pretty_table(data;
    column_labels=["Ambient", "Bundle E", "rk(E)", "χ(O_Z)", "h¹¹", "h²¹"],
    alignment=[:l, :l, :r, :r, :r, :r],
    fit_table_in_display_vertically=false,
    fit_table_in_display_horizontally=false,
  )
end

main()
