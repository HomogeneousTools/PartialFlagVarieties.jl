# ═══════════════════════════════════════════════════════════════════════════════
#  Benchmarks for PartialFlagVarieties.jl — with regression tracking
#
#  Usage:
#    julia --project=. benchmark/bench.jl               # run + save results
#    julia --project=. benchmark/bench.jl --compare      # compare vs last saved
#    julia --project=. benchmark/bench.jl --save-only    # save without comparing
#
#  Results are saved to benchmark/results/<timestamp>.json
#  The most recent saved result is also copied as benchmark/results/latest.json
# ═══════════════════════════════════════════════════════════════════════════════

using BenchmarkTools
using Printf
using Dates
using PartialFlagVarieties
using Lie

import Lie: borel_weil_bott

# ─── CLI ────────────────────────────────────────────────────────────────────

const COMPARE = "--compare" in ARGS
const SAVE_ONLY = "--save-only" in ARGS

# ─── Result storage ─────────────────────────────────────────────────────────

struct BenchResult
  name::String
  category::String
  min_time_ns::Float64
  median_time_ns::Float64
  mean_time_ns::Float64
  min_allocs::Int
  min_memory::Int      # bytes
  median_allocs::Int
  median_memory::Int   # bytes
end

const ALL_RESULTS = BenchResult[]

# ─── Helpers ─────────────────────────────────────────────────────────────────

function header(name)
  println("\n", "="^80)
  println("  ", name)
  println("="^80)
end

function report(name, b::BenchmarkTools.Trial; category="", extra="")
  t_min = minimum(b).time / 1e6
  t_med = median(b).time / 1e6
  t_mean = mean(b).time / 1e6
  a_min = minimum(b).allocs
  m_min = minimum(b).memory
  a_med = Int(round(median(b).allocs))
  m_med = Int(round(median(b).memory))

  @printf("  %-55s %9.3f ms  %6d allocs  %8s", name, t_min, a_min, fmt_bytes(m_min))
  extra != "" && print("  $extra")
  println()

  push!(
    ALL_RESULTS,
    BenchResult(name, category, minimum(b).time, median(b).time,
      mean(b).time, a_min, m_min, a_med, m_med),
  )
end

function fmt_bytes(b)
  b < 1024 && return @sprintf("%d B", b)
  b < 1024^2 && return @sprintf("%.1f KiB", b / 1024)
  b < 1024^3 && return @sprintf("%.1f MiB", b / 1024^2)
  return @sprintf("%.1f GiB", b / 1024^3)
end

"""Clear all PartialFlagVarieties + Lie caches for cold benchmarking."""
function clear_all_pfv_caches!()
  Lie.clear_all_caches!()
  PartialFlagVarieties.clear_caches!()
  return nothing
end

println("PartialFlagVarieties.jl Benchmarks")
println("="^80)

# ─── JIT warmup ─────────────────────────────────────────────────────────────
# Run a representative computation from each major code path to trigger
# compilation of all relevant methods *before* any benchmark timing begins.
# This prevents the first benchmark from absorbing JIT overhead that doesn't
# belong to it.

println("  Warming up JIT...")
let
  X = Gr(2, 4)
  T = tangent_bundle(X)
  Ω = cotangent_bundle(X)
  S = dual(universal_subbundle(X))
  E3 = exterior_power(T, 3)
  S2 = symmetric_power(S, 2)
  tp = tensor_product(S, universal_quotient_bundle(X))
  H = dimensions(E3)
  Z = zero_locus(reduce(direct_sum, [line_bundle(X, 1) for _ in 1:2]))
  hodge_numbers(Z)
  euler_characteristic(Z)
  cohomology_on_restriction(Z, structure_sheaf(X))
  hochschild_cohomology(X)
  beilinson_collection(X)
  is_strong_exceptional_sequence(beilinson_collection(X))
  filtered_tangent_bundle(X)
  schur_functor(X, (2, 1))
  is_calabi_yau_candidate(line_bundle(X, 4))
  # exceptional types
  Y = cayley_plane()
  betti_numbers(Y)
  tangent_bundle(Y)
  hodge_numbers(Y)
  # quadric path
  Q = quadric(3)
  spinor_bundle(Q)
  kapranov_collection(Q)
  # OGr path
  OGr(3, 6)
  SGr(2, 4)
end
clear_all_pfv_caches!()
println("  JIT warmup complete.\n")

# ═══════════════════════════════════════════════════════════════════════════════
#  1. Bundle algebra (cold — the real cost of decomposition)
# ═══════════════════════════════════════════════════════════════════════════════

header("1. Bundle algebra (cold)")

# Exterior powers of tangent bundle — the core tensor decomposition workload
# Each ∧ᵖT decomposes via multiexponents into IrrepLevi products calling Lie.jl

function _bench_ext_power_tangent(X, p)
  clear_all_pfv_caches!()
  exterior_power(tangent_bundle(X), p)
end

for (label, X, p) in [
  ("∧⁶T Gr(3,7)", Gr(3, 7), 6),
  ("∧⁸T Gr(3,7)", Gr(3, 7), 8),
  ("∧⁵T OGr(5,10)", OGr(5, 10), 5),
  ("∧⁴T adj G₂", adjoint_variety(TypeG2), 4),
  ("∧⁶T cayley_plane", cayley_plane(), 6),
  ("∧⁴T freudenthal", freudenthal_variety(), 4),
]
  _bench_ext_power_tangent(X, p)  # warmup
  b = @benchmark _bench_ext_power_tangent($X, $p) evals = 1 samples = 10
  report(label, b; category="bundle_algebra")
end

# Symmetric powers of dual universal subbundle (used in Küchle families)
function _bench_sym_power_subbundle(X, k)
  clear_all_pfv_caches!()
  symmetric_power(dual(universal_subbundle(X)), k)
end

for (label, X, k) in [
  ("Sym⁴S* Gr(2,7)", Gr(2, 7), 4),
  ("Sym³S* Gr(3,7)", Gr(3, 7), 3),
  ("Sym⁵S* Gr(2,7)", Gr(2, 7), 5),
]
  _bench_sym_power_subbundle(X, k)
  b = @benchmark _bench_sym_power_subbundle($X, $k) evals = 1 samples = 10
  report(label, b; category="bundle_algebra")
end

# Tensor product of non-trivial bundles (cold — triggers Lie.jl tensor_product)
function _bench_tensor_product(X)
  clear_all_pfv_caches!()
  S = dual(universal_subbundle(X))
  Q = universal_quotient_bundle(X)
  tensor_product(S, Q)
end

for (label, X) in [
  ("S*⊗Q Gr(2,7)", Gr(2, 7)),
  ("S*⊗Q Gr(3,7)", Gr(3, 7)),
]
  _bench_tensor_product(X)
  b = @benchmark _bench_tensor_product($X) evals = 1 samples = 10
  report(label, b; category="bundle_algebra")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  2. Cohomology — Borel–Weil–Bott on non-trivial bundles
# ═══════════════════════════════════════════════════════════════════════════════

header("2. Cohomology (cold)")

# Cohomology of exterior powers of tangent bundle (full BWB pipeline)
function _bench_cohom_ext_tangent(X, p)
  clear_all_pfv_caches!()
  dimensions(exterior_power(tangent_bundle(X), p))
end

for (label, X, p) in [
  ("H*(∧⁴T) Gr(2,5)", Gr(2, 5), 4),
  ("H*(∧⁶T) Gr(3,7)", Gr(3, 7), 6),
  ("H*(∧⁵T) OGr(5,10)", OGr(5, 10), 5),
  ("H*(∧⁵T) cayley_plane", cayley_plane(), 5),
]
  _bench_cohom_ext_tangent(X, p)
  b = @benchmark _bench_cohom_ext_tangent($X, $p) evals = 1 samples = 10
  report(label, b; category="cohomology")
end

# Bott vanishing search: ∧ᵖΩ(q) on adjoint varieties (from BottVanishing.jl)
function _bench_bott_vanishing_sweep(X, q_max)
  clear_all_pfv_caches!()
  Ω = cotangent_bundle(X)
  d = dimension(X)
  total = 0
  for p in 0:d
    Ωp = exterior_power(Ω, p)
    for q in 1:q_max
      Ωpq = twist(Ωp, 1, q)
      H = dimensions(Ωpq)
      total += sum(H[i] for i in 0:d)
    end
  end
  total
end

# B₃ adjoint (dim 9): sweep q=1..3 → 10×3 = 30 BWB calls
let X = adjoint_variety(TypeB{3})
  _bench_bott_vanishing_sweep(X, 3)
  b = @benchmark _bench_bott_vanishing_sweep($X, 3) evals = 1 samples = 3
  report("Bott sweep adj B₃ (q≤3)", b; category="cohomology")
end

# B₄ adjoint (dim 16): sweep q=1..2 → 17×2 = 34 BWB calls
let X = adjoint_variety(TypeB{4})
  _bench_bott_vanishing_sweep(X, 2)
  b = @benchmark _bench_bott_vanishing_sweep($X, 2) evals = 1 samples = 3
  report("Bott sweep adj B₄ (q≤2)", b; category="cohomology")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  3. Hodge numbers and Hochschild cohomology
# ═══════════════════════════════════════════════════════════════════════════════

header("3. Hodge / Hochschild (cold)")

# Hodge numbers on larger varieties (cold — forces full recomputation)
function _bench_hodge(X)
  clear_all_pfv_caches!()
  hodge_numbers(X)
end

for (label, X) in [
  ("cayley_plane", cayley_plane()),
  ("freudenthal_variety", freudenthal_variety()),
  ("OGr(5,10)", OGr(5, 10)),
]
  _bench_hodge(X)
  b = @benchmark _bench_hodge($X) evals = 1 samples = 5
  report("hodge_numbers($label)", b; category="hodge")
end

# Hochschild cohomology (polyvector parallelogram — many ext powers of T)
function _bench_hochschild(X)
  clear_all_pfv_caches!()
  hochschild_cohomology(X)
end

for (label, X) in [
  ("Gr(2,6)", Gr(2, 6)),
  ("Gr(3,6)", Gr(3, 6)),
  ("Gr(2,7)", Gr(2, 7)),
  ("cayley_plane", cayley_plane()),
]
  _bench_hochschild(X)
  b = @benchmark _bench_hochschild($X) evals = 1 samples = 5
  report("hochschild($label)", b; category="hodge")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  4. Exceptional collections (O(n²) pairwise checks)
# ═══════════════════════════════════════════════════════════════════════════════

header("4. Exceptional collections (cold)")

# Kapranov bundles on Gr(2,5): 10 bundles → 45 pair checks
function _bench_kapranov_exc(X)
  clear_all_pfv_caches!()
  coll = kapranov_bundles_grassmannian(X)
  is_strong_exceptional_sequence(coll)
end

let X = Gr(2, 5)
  _bench_kapranov_exc(X)
  b = @benchmark _bench_kapranov_exc($X) evals = 1 samples = 5
  report("strong_exc_seq(Kapranov Gr(2,5))", b; category="exceptional")
end

# Kapranov bundles on Gr(2,6): 15 bundles → 105 pair checks
let X = Gr(2, 6)
  _bench_kapranov_exc(X)
  b = @benchmark _bench_kapranov_exc($X) evals = 1 samples = 3
  report("strong_exc_seq(Kapranov Gr(2,6))", b; category="exceptional")
end

# Kapranov bundles on Gr(2,7): 21 bundles → 210 pair checks
let X = Gr(2, 7)
  _bench_kapranov_exc(X)
  b = @benchmark _bench_kapranov_exc($X) evals = 1 samples = 3
  report("strong_exc_seq(Kapranov Gr(2,7))", b; category="exceptional")
end

# Full exceptionality check on Gr(3,6): 20 bundles → 190 pair checks
function _bench_kapranov_full(X)
  clear_all_pfv_caches!()
  coll = kapranov_bundles_grassmannian(X)
  is_full_exceptional_sequence(coll, X)
end

let X = Gr(3, 6)
  _bench_kapranov_full(X)
  b = @benchmark _bench_kapranov_full($X) evals = 1 samples = 3
  report("full_exc_seq(Kapranov Gr(3,6))", b; category="exceptional")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  5. Zero loci — Hodge numbers via Koszul resolution
# ═══════════════════════════════════════════════════════════════════════════════

header("5. Zero loci — Hodge numbers")

# Quintic threefold: O(5) on ℙ⁴
function _bench_zl_hodge(E)
  clear_all_pfv_caches!()
  Z = zero_locus(E)
  hodge_numbers(Z)
end

let X = projective_space(4)
  E = line_bundle(X, 5)
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 10
  report("Hodge(quintic CY3 in ℙ⁴)", b; category="zero_loci")
end

# Fano threefold 1-15: O(1)³ on Gr(2,5)
let X = Gr(2, 5)
  E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:3])
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 5
  report("Hodge(O(1)³ on Gr(2,5))", b; category="zero_loci")
end

# Fano threefold 1-7: O(1)⁵ on Gr(2,6)  (rank 5, dim 8 ambient)
let X = Gr(2, 6)
  E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:5])
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 5
  report("Hodge(O(1)⁵ on Gr(2,6))", b; category="zero_loci")
end

# Fano lines on cubic fourfold: Sym³S* on Gr(2,6) — rank 4, dim 8
let X = Gr(2, 6)
  E = symmetric_power(dual(universal_subbundle(X)), 3)
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 5
  report("Hodge(Fano lines on cubic, Gr(2,6))", b; category="zero_loci")
end

# Küchle b9: (Sym²S*)² on Gr(2,7) — rank 6, dim 10
let X = Gr(2, 7)
  S = dual(universal_subbundle(X))
  E = direct_sum(symmetric_power(S, 2), symmetric_power(S, 2))
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 3
  report("Hodge(Küchle b9, Gr(2,7))", b; category="zero_loci")
end

# Linear section CY: 7·O(1) on Gr(2,7) — rank 7, dim 10
let X = Gr(2, 7)
  E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:7])
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 3
  report("Hodge(7·O(1) CY3 on Gr(2,7))", b; category="zero_loci")
end

# Fano threefold 1-10: (∧²U)³ on Gr(3,7) — rank 9, dim 12
let X = Gr(3, 7)
  U = universal_subbundle(X)
  E = reduce(direct_sum, [exterior_power(U, 2) for _ in 1:3])
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 3
  report("Hodge((∧²U)³ on Gr(3,7))", b; category="zero_loci")
end

# CY3 No.28 on Gr(3,8): (∧²S*)⁴ — rank 12, dim 15
let X = Gr(3, 8)
  S = dual(universal_subbundle(X))
  E = reduce(direct_sum, [exterior_power(S, 2) for _ in 1:4])
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 3
  report("Hodge((∧²S*)⁴ on Gr(3,8))", b; category="zero_loci")
end

# Hodge on OGr(5,10) CY3: 7·O(1)
let X = OGr(5, 10)
  E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:7])
  _bench_zl_hodge(E)
  b = @benchmark _bench_zl_hodge($E) evals = 1 samples = 3
  report("Hodge(7·O(1) CY3 on OGr(5,10))", b; category="zero_loci")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  6. Zero loci — Koszul cohomology on restriction
# ═══════════════════════════════════════════════════════════════════════════════

header("6. Zero loci — cohomology on restriction")

# Küchle-style: fano index computation via cohomology_on_restriction
function _bench_cohom_on_restriction(E, F)
  clear_all_pfv_caches!()
  Z = zero_locus(E)
  cohomology_on_restriction(Z, F)
end

# O(1)³ on Gr(2,5), restrict O(k) for k=1..3
let X = Gr(2, 5)
  E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:3])
  for k in [1, 2, 3]
    F = line_bundle(X, k)
    _bench_cohom_on_restriction(E, F)
    b = @benchmark _bench_cohom_on_restriction($E, $F) evals = 1 samples = 5
    report("H*(O($k))|_Z, O(1)³ Gr(2,5)", b; category="koszul")
  end
end

# Sym³S* on Gr(2,6), restrict O(k)
let X = Gr(2, 6)
  E = symmetric_power(dual(universal_subbundle(X)), 3)
  for k in [1, 2]
    F = line_bundle(X, k)
    _bench_cohom_on_restriction(E, F)
    b = @benchmark _bench_cohom_on_restriction($E, $F) evals = 1 samples = 3
    report("H*(O($k))|_Z, Sym³S* Gr(2,6)", b; category="koszul")
  end
end

# Küchle b9 on Gr(2,7), restrict O(1)
let X = Gr(2, 7)
  S = dual(universal_subbundle(X))
  E = direct_sum(symmetric_power(S, 2), symmetric_power(S, 2))
  F = line_bundle(X, 1)
  _bench_cohom_on_restriction(E, F)
  b = @benchmark _bench_cohom_on_restriction($E, $F) evals = 1 samples = 3
  report("H*(O(1))|_Z, Küchle b9 Gr(2,7)", b; category="koszul")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  7. Zero loci on exceptional varieties
# ═══════════════════════════════════════════════════════════════════════════════

header("7. Zero loci on exceptional varieties")

# CY3 on cayley plane: O(1)¹³ on OP² — heavy E₆ computation
function _bench_cayley_cy(E)
  clear_all_pfv_caches!()
  Z = zero_locus(E)
  euler_characteristic(Z)
end

let X = cayley_plane()
  E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:13])
  _bench_cayley_cy(E)
  b = @benchmark _bench_cayley_cy($E) evals = 1 samples = 3
  report("χ(CY3 in cayley_plane)", b; category="zero_loci_exceptional")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  8. End-to-end: Hochschild on G/P (from HochschildAffine.jl)
# ═══════════════════════════════════════════════════════════════════════════════

header("8. End-to-end: Hochschild on G/P")

# Full computation: for each p=0..dim, compute χ(∧ᵖT) via BWB
function _bench_hochschild_e2e(X)
  clear_all_pfv_caches!()
  T = tangent_bundle(X)
  d = dimension(X)
  total = BigInt(0)
  for p in 0:d
    Ep = exterior_power(T, p)
    H = dimensions(Ep)
    for i in 0:d
      total += H[i]
    end
  end
  total
end

for (label, X) in [
  ("Gr(2,5) dim=6", Gr(2, 5)),
  ("Gr(2,6) dim=8", Gr(2, 6)),
  ("Gr(3,6) dim=9", Gr(3, 6)),
  ("OGr(5,10) dim=10", OGr(5, 10)),
  ("cayley_plane dim=16", cayley_plane()),
]
  _bench_hochschild_e2e(X)
  b = @benchmark _bench_hochschild_e2e($X) evals = 1 samples = 3
  report("Hochschild e2e $label", b; category="end_to_end")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Save results
# ═══════════════════════════════════════════════════════════════════════════════

function save_results(results::Vector{BenchResult})
  results_dir = joinpath(@__DIR__, "results")
  mkpath(results_dir)

  timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
  filepath = joinpath(results_dir, "$timestamp.json")
  latest = joinpath(results_dir, "latest.json")

  open(filepath, "w") do io
    println(io, "{")
    println(io, "  \"timestamp\": \"$timestamp\",")
    println(io, "  \"julia_version\": \"$(VERSION)\",")
    println(io, "  \"results\": [")
    for (i, r) in enumerate(results)
      comma = i < length(results) ? "," : ""
      println(io, "    {")
      println(io, "      \"name\": $(repr(r.name)),")
      println(io, "      \"category\": $(repr(r.category)),")
      println(io, "      \"min_time_ns\": $(r.min_time_ns),")
      println(io, "      \"median_time_ns\": $(r.median_time_ns),")
      println(io, "      \"mean_time_ns\": $(r.mean_time_ns),")
      println(io, "      \"min_allocs\": $(r.min_allocs),")
      println(io, "      \"min_memory\": $(r.min_memory),")
      println(io, "      \"median_allocs\": $(r.median_allocs),")
      println(io, "      \"median_memory\": $(r.median_memory)")
      println(io, "    }$comma")
    end
    println(io, "  ]")
    println(io, "}")
  end

  isfile(latest) && rm(latest)
  cp(filepath, latest)

  println("\nResults saved to: $filepath")
  return filepath
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Compare against previous results
# ═══════════════════════════════════════════════════════════════════════════════

function load_results_json(filepath::String)
  content = read(filepath, String)
  results = Dict{String,Float64}()
  for m in
      eachmatch(r"\"name\":\s*\"([^\"]+)\"[^}]*\"min_time_ns\":\s*([0-9.e+\-]+)", content)
    results[m.captures[1]] = parse(Float64, m.captures[2])
  end
  return results
end

function compare_results_from_data(
  current::Vector{BenchResult}, baseline::Dict{String,Float64}, baseline_path::String
)
  println("\n", "="^80)
  println("  Regression comparison vs: $baseline_path")
  println("="^80)
  @printf("  %-55s %10s  %10s  %8s\n", "Benchmark", "Baseline", "Current", "Ratio")
  println("  ", "-"^90)

  regressions = 0
  improvements = 0

  for r in current
    if haskey(baseline, r.name)
      old_ns = baseline[r.name]
      new_ns = r.min_time_ns
      ratio = new_ns / old_ns

      marker = if ratio > 1.15
        regressions += 1
        "  REGRESSION"
      elseif ratio < 0.85
        improvements += 1
        "  IMPROVED"
      else
        ""
      end

      @printf("  %-55s %9.3f ms  %9.3f ms  %7.2fx%s\n",
        r.name, old_ns / 1e6, new_ns / 1e6, ratio, marker)
    else
      @printf("  %-55s %10s  %9.3f ms  %8s\n",
        r.name, "NEW", r.min_time_ns / 1e6, "-")
    end
  end

  println()
  println("  Summary: $improvements improved, $regressions regressions, ",
    "$(length(current) - improvements - regressions) unchanged")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

println("\n", "="^80)
println("  All benchmarks complete. $(length(ALL_RESULTS)) benchmarks recorded.")
println("="^80)

# Memory summary
total_allocs = sum(r.min_allocs for r in ALL_RESULTS)
total_memory = sum(r.min_memory for r in ALL_RESULTS)
println("  Total min allocations: $total_allocs")
println("  Total min memory: $(fmt_bytes(total_memory))")

# Per-category summary
categories = unique(r.category for r in ALL_RESULTS)
println("\n  Per-category timing (sum of min times):")
for cat in categories
  cat_results = filter(r -> r.category == cat, ALL_RESULTS)
  cat_time = sum(r.min_time_ns for r in cat_results) / 1e6
  cat_mem = sum(r.min_memory for r in cat_results)
  @printf("    %-30s %9.1f ms  %8s  (%d benchmarks)\n",
    cat, cat_time, fmt_bytes(cat_mem), length(cat_results))
end

# Load baseline BEFORE saving new results
baseline_path = if COMPARE
  lp = joinpath(@__DIR__, "results", "latest.json")
  isfile(lp) ? lp : nothing
else
  nothing
end
baseline_data = isnothing(baseline_path) ? nothing : load_results_json(baseline_path)

saved_path = save_results(ALL_RESULTS)

if COMPARE
  if !isnothing(baseline_data)
    compare_results_from_data(ALL_RESULTS, baseline_data, baseline_path)
  else
    println("No previous results to compare against.")
  end
end
