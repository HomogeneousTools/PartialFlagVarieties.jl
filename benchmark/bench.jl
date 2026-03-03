using BenchmarkTools
using PartialFlagVarieties
using Lie

println("PartialFlagVarieties.jl Benchmarks")
println("=" ^ 60)

# ─── Clear caches for fair benchmarking ──────────────────────────────────────
Lie.clear_all_caches!()

# ═══════════════════════════════════════════════════════════════════════════════
#  Type-level computation benchmarks (@generated functions)
# ═══════════════════════════════════════════════════════════════════════════════

println("\n── @generated function benchmarks (first-call + cached) ──\n")

suite_generated = BenchmarkGroup()

# Dimension
suite_generated["dimension/A4_P2"] = @benchmarkable dimension($MarkedDynkinType{TypeA{4},(2,)})
suite_generated["dimension/E6_P1"] = @benchmarkable dimension($MarkedDynkinType{TypeE{6},(1,)})
suite_generated["dimension/D5_P5"] = @benchmarkable dimension($MarkedDynkinType{TypeD{5},(5,)})

# Levi type
suite_generated["levi_type/A4_P2"] = @benchmarkable levi_type($MarkedDynkinType{TypeA{4},(2,)})
suite_generated["levi_type/E6_P1"] = @benchmarkable levi_type($MarkedDynkinType{TypeE{6},(1,)})

# Euler characteristic
suite_generated["euler_char/A4_P2"] = @benchmarkable euler_characteristic($MarkedDynkinType{TypeA{4},(2,)})
suite_generated["euler_char/E6_P1"] = @benchmarkable euler_characteristic($MarkedDynkinType{TypeE{6},(1,)})

# Betti numbers
suite_generated["betti/A4_P2"] = @benchmarkable betti_numbers($MarkedDynkinType{TypeA{4},(2,)})
suite_generated["betti/E6_P1"] = @benchmarkable betti_numbers($MarkedDynkinType{TypeE{6},(1,)})

# Special matrix
suite_generated["special_matrix/A3_P2"] = @benchmarkable special_matrix($MarkedDynkinType{TypeA{3},(2,)})
suite_generated["special_matrix/E6_P1"] = @benchmarkable special_matrix($MarkedDynkinType{TypeE{6},(1,)})

results_gen = run(suite_generated, seconds=2)
display(results_gen)

# ═══════════════════════════════════════════════════════════════════════════════
#  IrrepLevi benchmarks
# ═══════════════════════════════════════════════════════════════════════════════

println("\n── IrrepLevi benchmarks ──\n")

suite_levi = BenchmarkGroup()

let MDT = MarkedDynkinType{TypeA{4},(2,)}
  ω₁ = fundamental_weight(TypeA{4}, 1)
  ω₂ = fundamental_weight(TypeA{4}, 2)

  suite_levi["construct/A4_P2_ω₁"] = @benchmarkable IrrepLevi($MDT, $ω₁)
  suite_levi["construct/A4_P2_ω₂"] = @benchmarkable IrrepLevi($MDT, $ω₂)

  rep = IrrepLevi(MDT, ω₁)
  suite_levi["to_ambient/A4_P2"] = @benchmarkable to_ambient_weight($MDT, $rep)
  suite_levi["dual/A4_P2"] = @benchmarkable dual($rep)
  suite_levi["fiber_dim/A4_P2"] = @benchmarkable fiber_dimension($rep)
end

results_levi = run(suite_levi, seconds=2)
display(results_levi)

# ═══════════════════════════════════════════════════════════════════════════════
#  Bundle benchmarks
# ═══════════════════════════════════════════════════════════════════════════════

println("\n── Bundle benchmarks ──\n")

suite_bundle = BenchmarkGroup()

let MDT = MarkedDynkinType{TypeA{3},(2,)}
  suite_bundle["tangent/Gr24"] = @benchmarkable tangent_bundle($MDT)
  suite_bundle["cotangent/Gr24"] = @benchmarkable cotangent_bundle($MDT)
  suite_bundle["structure_sheaf/Gr24"] = @benchmarkable structure_sheaf($MDT)

  T = tangent_bundle(MDT)
  suite_bundle["rank/tangent_Gr24"] = @benchmarkable rank_bundle($T)
  suite_bundle["exterior_power_2/Gr24"] = @benchmarkable exterior_power($T, 2)
  suite_bundle["det/tangent_Gr24"] = @benchmarkable det_bundle($T)
end

let MDT = MarkedDynkinType{TypeA{4},(1,)}
  suite_bundle["tangent/P4"] = @benchmarkable tangent_bundle($MDT)
  T = tangent_bundle(MDT)
  suite_bundle["symmetric_power_2/P4"] = @benchmarkable symmetric_power($T, 2)
end

results_bundle = run(suite_bundle, seconds=2)
display(results_bundle)

# ═══════════════════════════════════════════════════════════════════════════════
#  Cohomology benchmarks
# ═══════════════════════════════════════════════════════════════════════════════

println("\n── Cohomology benchmarks ──\n")

suite_cohom = BenchmarkGroup()

let MDT = MarkedDynkinType{TypeA{4},(1,)}
  O = structure_sheaf(MDT)
  L = line_bundle(MDT, 1)

  suite_cohom["H*(P4, O)"] = @benchmarkable cohomology($MDT, $O)
  suite_cohom["H*(P4, O(1))"] = @benchmarkable cohomology($MDT, $L)
  suite_cohom["dim H*(P4, O)"] = @benchmarkable dimensions($MDT, $O)
end

let MDT = MarkedDynkinType{TypeA{3},(2,)}
  O = structure_sheaf(MDT)
  suite_cohom["H*(Gr24, O)"] = @benchmarkable cohomology($MDT, $O)
  suite_cohom["dim H*(Gr24, O)"] = @benchmarkable dimensions($MDT, $O)
end

results_cohom = run(suite_cohom, seconds=3)
display(results_cohom)

println("\nBenchmarks complete.")
