#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════════════════════
#  ExceptionalCollections.jl — examples for exceptional collection infrastructure
#
#  Run with:
#    julia --project=. examples/ExceptionalCollections.jl
#
#  Demonstrates:
#    1. Exceptionality predicates on projective space
#    2. Beilinson's exceptional collection on ℙ⁴
#    3. Dual Beilinson collection on ℙ⁴
#    4. Kapranov's collection on quadrics Q³, Q⁴, Q⁵
#    5. Schur functor ranks on Gr(k,n)
#    6. Kapranov–Orlov collection on Grassmannians Gr(2,4), Gr(2,5), Gr(3,6)
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables

# ─── 1. Exceptionality predicates on ℙ³ ────────────────────────────────────

println("=" ^ 70)
println("  1. Exceptionality predicates on ℙ³")
println("=" ^ 70)

X = projective_space(3)
bundles = [
  ("O",     structure_sheaf(X)),
  ("O(1)",  line_bundle(X, 1)),
  ("O(2)",  line_bundle(X, 2)),
  ("O(-1)", line_bundle(X, -1)),
  ("T",     tangent_bundle(X)),
  ("Ω",     cotangent_bundle(X)),
]

data = [Any[name, is_exceptional(E)] for (name, E) in bundles]
data_matrix = permutedims(hcat(data...), (2, 1))
pretty_table(
  data_matrix;
  column_labels = ["Bundle", "is_exceptional"],
  alignment = [:l, :c],
)

# ─── 2. Beilinson collection on ℙ⁴ ─────────────────────────────────────────

println("=" ^ 70)
println("  2. Beilinson collection on ℙ⁴")
println("=" ^ 70)

X4 = projective_space(4)
Es = beilinson_collection(X4)
println("Collection: ⟨O, O(1), O(2), O(3), O(4)⟩")
println("  Length          : ", length(Es))
println("  χ(ℙ⁴)           : ", euler_characteristic(X4))
println("  full exc. seq.  : ", is_full_exceptional_sequence(Es, X4))
println("  strong exc. seq.: ", is_strong_exceptional_sequence(Es))
println()

# ─── 3. Dual Beilinson collection on ℙ⁴ ────────────────────────────────────

println("=" ^ 70)
println("  3. Dual Beilinson collection on ℙ⁴")
println("=" ^ 70)

Ed = beilinson_collection_dual(X4)
bundle_names = ["Ω⁴(4)", "Ω³(3)", "Ω²(2)", "Ω¹(1)", "O"]
println("Collection: ⟨Ω⁴(4), Ω³(3), Ω²(2), Ω(1), O⟩")
println("  Ranks           : ", [rank_bundle(E) for E in Ed])
println("  full exc. seq.  : ", is_full_exceptional_sequence(Ed, X4))
println("  strong exc. seq.: ", is_strong_exceptional_sequence(Ed))
println()

# ─── 4. Kapranov collection on quadrics ─────────────────────────────────────

println("=" ^ 70)
println("  4. Kapranov collection on quadrics")
println("=" ^ 70)

quadric_data = Matrix{Any}(undef, 3, 5)
for (i, n) in enumerate([3, 4, 5])
  local Q, Es, chi, full, strong
  Q = quadric(n)
  Es = kapranov_collection(Q)
  chi = euler_characteristic(Q)
  full = is_full_exceptional_sequence(Es, Q)
  strong = is_strong_exceptional_sequence(Es)
  quadric_data[i, :] = ["Q^$n", chi, length(Es), full, strong]
end
pretty_table(
  quadric_data;
  column_labels = ["Variety", "χ", "length", "full exc.", "strong exc."],
  alignment = [:l, :c, :c, :c, :c],
)

# ─── 5. Schur functor ranks on Gr(k,n) ──────────────────────────────────────

println("=" ^ 70)
println("  5. Schur functor Σ^α(U^∨) ranks on Gr(2,5)")
println("=" ^ 70)

X25 = Gr(2, 5)
partitions25 = [[a, b] for a in 0:3 for b in 0:a]
schur_data = [begin
  E = schur_functor(X25, α)
  [string(α), rank_bundle(E), is_exceptional(E)]
end for α in partitions25]
schur_matrix = permutedims(hcat(schur_data...), (2, 1))
pretty_table(
  schur_matrix;
  column_labels = ["partition α", "rank(Σ^α U^∨)", "is_exceptional"],
  alignment = [:l, :c, :c],
)

# ─── 6. Kapranov–Orlov collection on Grassmannians ──────────────────────────

println("=" ^ 70)
println("  6. Kapranov–Orlov collection on Grassmannians")
println("=" ^ 70)

grass_data = Matrix{Any}(undef, 4, 5)
for (i, (k, n)) in enumerate([(2, 4), (2, 5), (2, 6), (3, 6)])
  local X, Es, chi, full, strong
  X = Gr(k, n)
  Es = kapranov_bundles_grassmannian(X)
  chi = euler_characteristic(X)
  full = is_full_exceptional_sequence(Es, X)
  strong = is_strong_exceptional_sequence(Es)
  grass_data[i, :] = ["Gr($k,$n)", chi, length(Es), full, strong]
end
pretty_table(
  grass_data;
  column_labels = ["Variety", "χ", "length", "full exc.", "strong exc."],
  alignment = [:l, :c, :c, :c, :c],
)
