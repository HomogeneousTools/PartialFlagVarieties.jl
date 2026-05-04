#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════════════════════
#  HyperellipticLinesExceptional.jl — Exceptional collections on F₁(Q₁ ∩ Q₂)
#
#  Verifies the conjectured exceptional collection on the Fano scheme of
#  lines F₁(Q₁ ∩ Q₂) on the intersection of two quadrics in P^{2g+1},
#  realized as the zero locus of E = (Sym²U^∨)^⊕2 on Gr(2, 2g+2).
#
#  The conjecture predicts (g-1)(2g-5) exceptional objects with a Lefschetz
#  decomposition whose initial block is
#    ⟨O, U^∨, Sym²U^∨, …, Sym^{g-2}U^∨⟩
#  and support partition (g-1,…,g-1; g-2,…,g-2) with (g-3) copies of (g-1)
#  and (g-1) copies of (g-2).
#
#  Also verifies the stacky variant on Gr(2, 2g+1) with (g-2)(2g-5) objects
#  and a rectangular Lefschetz decomposition.
#
#  Reference:
#    Exceptional objects on F₁(Q₁ ∩ Q₂) — hyperelliptic-lines-exceptional.tex
#
#  Run with:
#    julia --project=. examples/HyperellipticLinesExceptional.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables

# ═══════════════════════════════════════════════════════════════════════════════
#  Build the conjectured exceptional collection for genus g
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hyperelliptic_collection(g)

Build the conjectured exceptional collection on F₁(Q₁ ∩ Q₂) ⊂ Gr(2, 2g+2).

Returns `(G, E, Z, L)` where:
- `G = Gr(2, 2g+2)` is the ambient Grassmannian,
- `E = (Sym²U^∨)^⊕2` is the defining bundle,
- `Z = Z(E)` is the zero locus,
- `L` is the list of (g-1)(2g-5) exceptional objects (bundles on G).

The collection has a Lefschetz structure:
- For k = 0, …, g-4: full block {Sym^i(U^∨)(k) : i = 0, …, g-2}
- For k = g-3, …, 2g-5: short block {Sym^i(U^∨)(k) : i = 0, …, g-3}
"""
function hyperelliptic_collection(g::Int)
  G = Gr(2, 2 * g + 2)
  U = universal_subbundle(G)
  E = 2 * symmetric_power(dual(U), 2)
  Z = zero_locus(E)

  L = CompletelyReducibleBundle[]
  for k in 0:(2 * g - 5)
    for i in 0:(g - 2)
      # Skip the top symmetric power in the short blocks
      if k >= g - 3 && i == g - 2
        continue
      end
      F = if i == 0
        twist(structure_sheaf(G), 1, k)
      else
        twist(symmetric_power(dual(U), i), 1, k)
      end
      push!(L, F)
    end
  end

  (G, E, Z, L)
end

"""
    hyperelliptic_collection_stacky(g)

Build the conjectured rectangular Lefschetz exceptional collection on the
stacky variant, realized as Z((Sym²U^∨)^⊕2) on Gr(2, 2g+1).

The collection has (g-2)(2g-5) objects:
  {Sym^i(U^∨)(k) : i = 0, …, g-3,  k = 0, …, 2g-6}
"""
function hyperelliptic_collection_stacky(g::Int)
  G = Gr(2, 2 * g + 1)
  U = universal_subbundle(G)
  E = 2 * symmetric_power(dual(U), 2)
  Z = zero_locus(E)

  L = CompletelyReducibleBundle[]
  for k in 0:(2 * g - 6)
    for i in 0:(g - 3)
      F = if i == 0
        twist(structure_sheaf(G), 1, k)
      else
        twist(symmetric_power(dual(U), i), 1, k)
      end
      push!(L, F)
    end
  end

  (G, E, Z, L)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Check the conjecture
# ═══════════════════════════════════════════════════════════════════════════════

function check_conjecture(g::Int)
  println("─" ^ 70)
  println("  g = $g:  F₁(Q₁ ∩ Q₂) = Z((Sym²U^∨)^⊕2) ⊂ Gr(2, $(2g+2))")
  println("─" ^ 70)

  (G, E, Z, L) = hyperelliptic_collection(g)

  expected = (g - 1) * (2 * g - 5)
  println("  dim(Gr)   = ", dimension(G))
  println("  rank(E)   = ", rank_bundle(E))
  println("  dim(Z)    = ", dimension(Z))
  println("  expected  = (g-1)(2g-5) = $expected objects")
  println("  actual    = ", length(L), " objects")
  @assert length(L) == expected "Expected $expected objects, got $(length(L))"

  print("  checking exceptional sequence... ")
  exc = is_exceptional_sequence(L, Z)
  println(exc ? "✓" : "✗")

  print("  checking strong exceptional sequence... ")
  strong = is_strong_exceptional_sequence(L, Z)
  println(strong ? "✓" : "✗")

  println()
  (exc, strong)
end

function check_conjecture_stacky(g::Int)
  println("─" ^ 70)
  println("  g = $g (stacky):  Z((Sym²U^∨)^⊕2) ⊂ Gr(2, $(2g+1))")
  println("─" ^ 70)

  (G, E, Z, L) = hyperelliptic_collection_stacky(g)

  expected = (g - 2) * (2 * g - 5)
  println("  dim(Gr)   = ", dimension(G))
  println("  rank(E)   = ", rank_bundle(E))
  println("  dim(Z)    = ", dimension(Z))
  println("  expected  = (g-2)(2g-5) = $expected objects")
  println("  actual    = ", length(L), " objects")
  @assert length(L) == expected "Expected $expected objects, got $(length(L))"

  print("  checking exceptional sequence... ")
  exc = is_exceptional_sequence(L, Z)
  println(exc ? "✓" : "✗")

  print("  checking strong exceptional sequence... ")
  strong = is_strong_exceptional_sequence(L, Z)
  println(strong ? "✓" : "✗")

  println()
  (exc, strong)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

println("=" ^ 70)
println("  Exceptional collections on Fano schemes of lines F₁(Q₁ ∩ Q₂)")
println("=" ^ 70)
println()

results = Any[]

for g in 3:8
  (exc, strong) = check_conjecture(g)
  push!(results, ["g=$g", "Gr(2,$(2g+2))", (g - 1) * (2g - 5), exc, strong])
end

for g in 3:8
  (exc, strong) = check_conjecture_stacky(g)
  push!(results, ["g=$g stacky", "Gr(2,$(2g+1))", (g - 2) * (2g - 5), exc, strong])
end

# ─── Summary table ──────────────────────────────────────────────────────────

println("=" ^ 70)
println("  Summary")
println("=" ^ 70)

data = permutedims(hcat(results...), (2, 1))
pretty_table(
  data;
  column_labels=["Case", "Ambient", "# objects", "exceptional", "strong"],
  alignment=[:l, :c, :c, :c, :c],
)
