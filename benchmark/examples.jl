# ═══════════════════════════════════════════════════════════════════════════════
#  Example-level benchmarks for PartialFlagVarieties.jl
#
#  Runs each example script as a subprocess, captures wall-clock time and exit
#  codes. Repeats each example N times for statistical stability.
#
#  Usage:
#    julia --project=. benchmark/examples.jl                # run + save baseline
#    julia --project=. benchmark/examples.jl --compare       # compare vs baseline
#    julia --project=. benchmark/examples.jl --runs 5        # 5 repetitions
#    julia --project=. benchmark/examples.jl --tag pre-lru   # custom tag in filename
#
#  Results are saved to benchmark/results/examples_<tag>_<timestamp>.json
# ═══════════════════════════════════════════════════════════════════════════════

using Dates
using Printf

# ─── CLI parsing ──────────────────────────────────────────────────────────────

const COMPARE = "--compare" in ARGS
const N_RUNS = let
  idx = findfirst(==("--runs"), ARGS)
  idx !== nothing && idx < length(ARGS) ? parse(Int, ARGS[idx + 1]) : 3
end
const TAG = let
  idx = findfirst(==("--tag"), ARGS)
  idx !== nothing && idx < length(ARGS) ? ARGS[idx + 1] : ""
end

# ─── Example definitions ─────────────────────────────────────────────────────

const PROJECT_ROOT = dirname(@__DIR__)
const EXAMPLES_DIR = joinpath(PROJECT_ROOT, "examples")
const RESULTS_DIR = joinpath(@__DIR__, "results")

# Examples to benchmark: (name, filename, extra_args)
# Excluded: HochschildAffine.jl (too slow for repeated runs)
const EXAMPLES = [
  ("Hyperkaehler", "Hyperkaehler.jl", String[]),
  ("ExceptionalCollections", "ExceptionalCollections.jl", String[]),
  ("BottVanishing", "BottVanishing.jl", String[]),
  ("HyperellipticLinesExceptional", "HyperellipticLinesExceptional.jl", String[]),
  ("LinearSections", "LinearSections.jl", String[]),
  ("FanoThreefolds", "FanoThreefolds.jl", String[]),
  ("Kuechle", "Kuechle.jl", String[]),
  ("CICY3-1606.04076", "CICY3-1606.04076.jl", String[]),
  ("CICY3-1607.07821", "CICY3-1607.07821.jl", String[]),
  ("FanoFourfolds", "FanoFourfolds.jl", String[]),
  ("FanoFourfoldsK3", "FanoFourfoldsK3.jl", String[]),
  ("FanoFourfoldsQuiverZeroLoci", "FanoFourfoldsQuiverZeroLoci.jl", String[]),
]

# ─── Helpers ──────────────────────────────────────────────────────────────────

function fmt_time(seconds)
  if seconds < 60
    @sprintf("%.1fs", seconds)
  elseif seconds < 3600
    m = div(seconds, 60)
    s = seconds - m * 60
    @sprintf("%dm %.1fs", m, s)
  else
    h = div(seconds, 3600)
    rest = seconds - h * 3600
    m = div(rest, 60)
    @sprintf("%dh %dm", h, m)
  end
end

function run_example(name::String, filename::String, extra_args::Vector{String})
  script = joinpath(EXAMPLES_DIR, filename)
  if !isfile(script)
    @warn "Example not found: $script"
    return (exit_code=-1, elapsed=NaN)
  end

  cmd = `$(Base.julia_cmd()) --project=$PROJECT_ROOT $script $extra_args`

  t0 = time()
  proc = run(pipeline(cmd; stdout=devnull, stderr=devnull); wait=true)
  elapsed = time() - t0

  (exit_code=proc.exitcode, elapsed=elapsed)
end

# ─── Main benchmark loop ─────────────────────────────────────────────────────

function run_benchmarks(; n_runs::Int=N_RUNS)
  println("PartialFlagVarieties.jl — Example Benchmarks")
  println("="^80)
  println("  Examples:   $(length(EXAMPLES))")
  println("  Runs each:  $n_runs")
  println("  Julia:      $(VERSION)")
  println("="^80)
  println()

  results = []

  for (i, (name, filename, extra_args)) in enumerate(EXAMPLES)
    @printf("  [%2d/%d] %-40s", i, length(EXAMPLES), name)
    flush(stdout)

    times = Float64[]
    exit_codes = Int[]

    for run_i in 1:n_runs
      r = run_example(name, filename, extra_args)
      push!(exit_codes, r.exit_code)
      push!(times, r.elapsed)

      if run_i < n_runs
        print(".")
        flush(stdout)
      end
    end

    med = sort(times)[div(length(times), 2) + 1]
    mn = minimum(times)
    status = all(==(0), exit_codes) ? "OK" : "FAIL"

    @printf("  %8s  min=%8s  med=%8s  [%s]\n",
      status, fmt_time(mn), fmt_time(med),
      join([fmt_time(t) for t in times], ", "))

    push!(
      results,
      Dict(
        "name" => name,
        "filename" => filename,
        "times_s" => times,
        "min_s" => mn,
        "median_s" => med,
        "exit_codes" => exit_codes,
        "status" => status,
      ),
    )
  end

  results
end

# ─── Save / Load results ─────────────────────────────────────────────────────

function save_results(results; tag::String=TAG)
  mkpath(RESULTS_DIR)

  timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
  tag_part = isempty(tag) ? "" : "_$(tag)"
  filename = "examples$(tag_part)_$(timestamp).json"
  filepath = joinpath(RESULTS_DIR, filename)
  latest = joinpath(RESULTS_DIR, "examples_latest.json")

  open(filepath, "w") do io
    println(io, "{")
    println(io, "  \"timestamp\": \"$timestamp\",")
    println(io, "  \"julia_version\": \"$(VERSION)\",")
    println(io, "  \"tag\": $(repr(tag)),")
    println(io, "  \"n_runs\": $N_RUNS,")
    println(io, "  \"results\": [")
    for (i, r) in enumerate(results)
      comma = i < length(results) ? "," : ""
      println(io, "    {")
      println(io, "      \"name\": $(repr(r["name"])),")
      println(io, "      \"filename\": $(repr(r["filename"])),")
      println(
        io,
        "      \"times_s\": [$(join([@sprintf("%.3f", t) for t in r["times_s"]], ", "))],",
      )
      println(io, "      \"min_s\": $(@sprintf("%.3f", r["min_s"])),")
      println(io, "      \"median_s\": $(@sprintf("%.3f", r["median_s"])),")
      println(io, "      \"exit_codes\": [$(join(r["exit_codes"], ", "))],")
      println(io, "      \"status\": $(repr(r["status"]))")
      println(io, "    }$comma")
    end
    println(io, "  ]")
    println(io, "}")
  end

  isfile(latest) && rm(latest)
  cp(filepath, latest)

  println("\nResults saved to: $filepath")
  filepath
end

function load_results_json(filepath::String)
  content = read(filepath, String)
  results = Dict{String,Float64}()
  for m in eachmatch(r"\"name\":\s*\"([^\"]+)\"[^}]*\"min_s\":\s*([0-9.]+)", content)
    results[m.captures[1]] = parse(Float64, m.captures[2])
  end
  results
end

function compare_results(current, baseline_path::String)
  baseline = load_results_json(baseline_path)

  println("\n", "="^80)
  println("  Comparison vs: $(basename(baseline_path))")
  println("="^80)
  @printf("  %-40s %10s  %10s  %8s\n", "Example", "Baseline", "Current", "Ratio")
  println("  ", "-"^75)

  regressions = 0
  improvements = 0

  for r in current
    name = r["name"]
    new_t = r["min_s"]

    if haskey(baseline, name)
      old_t = baseline[name]
      ratio = new_t / old_t

      marker = if ratio > 1.15
        regressions += 1
        "  REGRESSION"
      elseif ratio < 0.85
        improvements += 1
        "  IMPROVED"
      else
        ""
      end

      @printf("  %-40s %9s  %9s  %7.2fx%s\n",
        name, fmt_time(old_t), fmt_time(new_t), ratio, marker)
    else
      @printf("  %-40s %10s  %9s  %8s\n",
        name, "NEW", fmt_time(new_t), "-")
    end
  end

  println()
  println("  Summary: $improvements improved, $regressions regressions, ",
    "$(length(current) - improvements - regressions) unchanged")
end

# ─── Main ─────────────────────────────────────────────────────────────────────

# Load baseline before running (if comparing)
baseline_path = if COMPARE
  lp = joinpath(RESULTS_DIR, "examples_latest.json")
  isfile(lp) ? lp : nothing
else
  nothing
end

results = run_benchmarks()
saved_path = save_results(results)

if COMPARE
  if !isnothing(baseline_path)
    compare_results(results, baseline_path)
  else
    println("No previous results to compare against.")
  end
end
