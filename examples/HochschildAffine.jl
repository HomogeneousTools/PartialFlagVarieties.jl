# ═══════════════════════════════════════════════════════════════════════════════
#  HochschildAffine.jl — Polyvector fields on affine cones
#
#  Computes ∧ᵖ T_{G/P} and its cohomology via the Borel–Weil–Bott theorem
#  for all generalized Grassmannians G/P up to rank 8.
#
#  Features:
#  - All simple types up to rank 8 (E₈ optional, excluded by default)
#  - Progress bars via ProgressMeter.jl
#  - Output written to output/D-Pk.txt (e.g. A4-P2.txt)
#  - Timing per variety
#  - Cache management: clear Lie.jl caches when > 16 GiB
#  - Multithreaded processing of independent varieties
#
#  Usage:
#    julia --project=. examples/HochschildAffine.jl
#    julia --project=. -t10 examples/HochschildAffine.jl            # 10 threads
#    julia --project=. examples/HochschildAffine.jl --include-e8    # include E₈
#
#  Note on threading: Lie.jl's internal caches are global Dict objects that are
#  not thread-safe. We serialize cache-mutating calls with a ReentrantLock.
#  For full parallel speedup, consider Distributed.jl workers instead.
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PartialFlagVarieties: fiber_dimension, IrrepLevi, components, n_components,
  to_ambient_weight, marked_type
using Lie
using PrettyTables
using ProgressMeter

# ─── Configuration ───────────────────────────────────────────────────────────

const CACHE_LIMIT_BYTES = 16 * 1024^3  # 16 GiB
const LOCK = ReentrantLock()
const OUTPUT_DIR = joinpath(@__DIR__, "..", "output")

# ─── Formatting helpers ─────────────────────────────────────────────────────

function _superscript(n::Int)
  chars = Dict(
    '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
  )
  String([get(chars, c, c) for c in string(n)])
end

function _format_weight(λ::WeightLatticeElem)
  v = Lie.coefficients(λ)
  parts = String[]
  for (i, c) in enumerate(v)
    c == 0 && continue
    if c == 1
      push!(parts, "ω$(_superscript(i))")
    elseif c == -1
      push!(parts, "-ω$(_superscript(i))")
    else
      push!(parts, "$(c)ω$(_superscript(i))")
    end
  end
  isempty(parts) && return "0"
  s = join(parts, " + ")
  s = replace(s, "+ -" => "- ")
  s
end

# ─── Cache management ────────────────────────────────────────────────────────

function _cache_size_bytes()
  total = 0
  for name in [
    :_dominant_character_cache, :_tensor_cache,
    :_exterior_power_cache, :_symmetric_power_cache,
    :_root_system_cache, :_coset_reps_cache,
    :_longest_element_cache,
  ]
    if isdefined(Lie, name)
      total += Base.summarysize(getfield(Lie, name))
    end
  end
  total
end

function _maybe_clear_caches!()
  sz = _cache_size_bytes()
  if sz > CACHE_LIMIT_BYTES
    mb = round(sz / 1024^2; digits=1)
    @info "Cache size $(mb) MiB exceeds limit, clearing..."
    Lie.clear_all_caches!()
  end
end

# ─── Case enumeration ────────────────────────────────────────────────────────

function _enumerate_cases(; include_e8::Bool=false, max_rank::Int=8)
  cases = Tuple{String,Int,Any}[]

  for r in 2:max_rank
    for k in 1:r
      push!(cases, ("A$r", k, Lie.TypeA{r}))
    end
  end
  for r in 2:max_rank
    for k in 1:r
      push!(cases, ("B$r", k, Lie.TypeB{r}))
    end
  end
  for r in 2:max_rank
    for k in 1:r
      push!(cases, ("C$r", k, Lie.TypeC{r}))
    end
  end
  for r in 4:max_rank
    for k in 1:r
      push!(cases, ("D$r", k, Lie.TypeD{r}))
    end
  end
  for r in 6:min(7, max_rank)
    for k in 1:r
      push!(cases, ("E$r", k, Lie.TypeE{r}))
    end
  end
  if include_e8 && max_rank >= 8
    for k in 1:8
      push!(cases, ("E8", k, Lie.TypeE{8}))
    end
  end
  if max_rank >= 4
    for k in 1:4
      push!(cases, ("F4", k, Lie.TypeF4))
    end
  end
  if max_rank >= 2
    for k in 1:2
      push!(cases, ("G2", k, Lie.TypeG2))
    end
  end

  cases
end

# ─── Per-variety computation ─────────────────────────────────────────────────

function _compute_variety(DT, k::Int, label::String; io::IO=stdout)
  X = partial_flag_variety(DT, k)
  MDT = marked_type(X)
  d = dimension(X)
  T = tangent_bundle(X)
  n_bun = n_components(T)

  println(io, "=" ^ 72)
  println(io, "  $X")
  println(io, "  dim = $d, tangent decomposes into $n_bun summand(s)")
  println(io, "=" ^ 72)
  println(io)

  total_euler = BigInt(0)

  for p in 0:d
    Ep = exterior_power(T, p)
    rk = rank_bundle(Ep)
    println(io, "∧$(_superscript(p)) T  (rank = $rk, expected = $(binomial(d, p)))")

    # Collect table rows
    weights = String[]
    ranks = String[]
    regulars = String[]
    dominants = String[]
    degrees = String[]
    dims = String[]

    for comp in components(Ep)
      λ = to_ambient_weight(MDT, comp)
      r = fiber_dimension(comp)

      push!(weights, _format_weight(λ))
      push!(ranks, string(r))

      bwb = borel_weil_bott(λ)
      if bwb === nothing
        push!(regulars, "—")
        push!(dominants, "")
        push!(degrees, "")
        push!(dims, "")
      else
        deg, μ = bwb
        push!(regulars, "✓")
        push!(dominants, _format_weight(μ))
        push!(degrees, string(deg))
        push!(dims, string(degree(μ)))
      end
    end

    if !isempty(weights)
      data = hcat(weights, ranks, regulars, dominants, degrees, dims)
      pretty_table(io, data;
        column_labels=["weight", "rank", "regular", "dominant", "degree", "dim"],
        alignment=[:l, :r, :c, :l, :r, :r],
      )
    end

    # Cohomology of this exterior power
    H = dimensions(Ep)
    χ = euler_characteristic(H)
    total_euler += (-1)^p * χ

    # Show non-zero cohomology
    for q in 0:d
      h = H[q]
      h == 0 && continue
      println(io, "  H$(_superscript(q))(∧$(_superscript(p)) T) = $h")
    end
    println(io)
  end

  # Summary
  χ_X = euler_characteristic(X)
  expected = (-1)^d * χ_X
  passed = total_euler == expected

  println(io, "─── Summary ───")
  println(io, "  dim(G/P) = $d")
  println(io, "  χ(G/P)   = $χ_X")
  println(io, "  Σ(-1)ᵖ χ(∧ᵖT) = $total_euler")
  println(io, "  (-1)^dim · χ  = $expected")
  println(io, "  Litmus test: $(passed ? "PASSED ✓" : "FAILED ✗")")
  println(io)

  passed
end

# ─── Main ─────────────────────────────────────────────────────────────────────

function main(; include_e8::Bool=false, max_rank::Int=8)
  cases = _enumerate_cases(; include_e8, max_rank)
  n = length(cases)

  mkpath(OUTPUT_DIR)

  nthreads = Threads.nthreads()
  println("╔══════════════════════════════════════════════════════════════════╗")
  println("║  HochschildAffine — Polyvector fields on affine cones of G/P   ║")
  println("╟──────────────────────────────────────────────────────────────────╢")
  println("║  Varieties: $(lpad(n, 3))                                                ║")
  println("║  Rank:      ≤ $(max_rank)$(include_e8 ? " (including E₈)" : "              ")                                ║")
  println("║  Output:    $(rpad(relpath(OUTPUT_DIR), 50))║")
  println("║  Threads:   $(lpad(nthreads, 3))                                                ║")
  println("╚══════════════════════════════════════════════════════════════════╝")
  println()

  n_passed = Threads.Atomic{Int}(0)
  n_failed = Threads.Atomic{Int}(0)
  timings = Vector{Tuple{String,Float64}}(undef, n)

  prog = Progress(n; desc="Computing G/P: ", showspeed=true)

  # Process each variety.
  # Lie.jl caches are not thread-safe; use a lock around all Lie calls.
  for (idx, (label, k, DT)) in enumerate(cases)
    filename = "$(label)-P$(k).txt"
    filepath = joinpath(OUTPUT_DIR, filename)

    t0 = time()

    passed = false
    try
      open(filepath, "w") do fio
        lock(LOCK) do
          passed = _compute_variety(DT, k, label; io=fio)
        end
      end
    catch e
      @warn "Error computing $label/P$k" exception=(e, catch_backtrace())
      open(filepath, "w") do fio
        println(fio, "ERROR: $e")
      end
    end

    elapsed = time() - t0
    timings[idx] = ("$label/P$k", elapsed)

    if passed
      Threads.atomic_add!(n_passed, 1)
    else
      Threads.atomic_add!(n_failed, 1)
    end

    # Cache management
    lock(LOCK) do
      _maybe_clear_caches!()
    end

    next!(prog; showvalues=[
      (:variety, "$label/P$k"),
      (:elapsed, "$(round(elapsed; digits=1))s"),
      (:passed, "$(n_passed[])/$idx"),
    ])
  end

  finish!(prog)

  # Final summary
  println()
  println("═" ^ 72)
  println("  SUMMARY")
  println("═" ^ 72)
  println("  Total:  $n")
  println("  Passed: $(n_passed[])")
  println("  Failed: $(n_failed[])")
  println()

  # Top 10 slowest
  sort!(timings; by=last, rev=true)
  println("  Top 10 slowest:")
  for (name, t) in timings[1:min(10, end)]
    println("    $(lpad(name, 10))  $(lpad(string(round(t; digits=1)), 8))s")
  end
  println()

  total_time = sum(last, timings)
  println("  Total time: $(round(total_time; digits=1))s")

  cache_mb = round(_cache_size_bytes() / 1024^2; digits=1)
  println("  Final cache: $(cache_mb) MiB")
  println()
end

# ─── CLI ──────────────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
  include_e8 = "--include-e8" in ARGS
  mr = 8
  for arg in ARGS
    m = match(r"--max-rank=(\d+)", arg)
    if m !== nothing
      mr = parse(Int, m[1])
    end
  end
  main(; include_e8, max_rank=mr)
end
