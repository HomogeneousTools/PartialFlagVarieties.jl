# ═══════════════════════════════════════════════════════════════════════════════
#  FanoFourfolds.jl — Hodge numbers of Fano fourfolds from data.json
#
#  Reads the dataset of Fano fourfold zero loci in products of type A
#  flag varieties (from Fatighenti–Mongardi et al.), constructs each
#  ambient variety and equivariant bundle in PartialFlagVarieties.jl,
#  computes Hodge numbers via the Koszul complex, and compares with the
#  reference values.
#
#  Data format (Schubert2/Macaulay2 convention):
#   - ambient: list of factors [k₁,...,kₘ, n], each encoding Fl(k₁,...,kₘ; n)
#   - bundle: nested arrays encoding direct sums of box products of
#     Schur functor bundles
#
#  Usage: julia --project=. examples/FanoFourfolds.jl [data.json path]
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using Lie
using JSON
using PrettyTables
using Printf
using ProgressMeter

# =============================================================================
#  GL(n) weight → omega coordinates (general flag)
# =============================================================================

"""
Convert a Schubert2-convention GL(n) weight vector to omega coordinates.

For a flag variety Fl(k₁,...,kₘ; n), the weight vector `w` of length `n`
stores eigenvalues in Schubert2 block order (quotient first, smallest
tautological sub last).  This function reorders to standard epsilon order
(smallest sub first) and takes successive differences.

# Arguments
- `ks`: sorted step sizes `[k₁, k₂, ..., kₘ]` with `k₁ < k₂ < ... < kₘ < n`
- `n`:  ambient dimension
- `w`:  GL(n) weight vector of length `n`
"""
function gl_weight_to_omega_flag(ks::Vector{Int}, n::Int, w::Vector{Int})
  # Schubert2 block order (left to right):
  #   Q*:           w[1 : n-kₘ]           (size n-kₘ)
  #   (Uₘ/Uₘ₋₁)*:  w[n-kₘ+1 : n-kₘ₋₁]  (size kₘ-kₘ₋₁)
  #   ...
  #   U₁*:          w[n-k₁+1 : n]         (size k₁)
  #
  # Standard epsilon order: [U₁*, (U₂/U₁)*, ..., Q*]
  eps = Int[]
  append!(eps, w[n - ks[1] + 1:n])            # U₁*
  for i in 2:length(ks)
    append!(eps, w[n - ks[i] + 1:n - ks[i - 1]])  # (Uᵢ/Uᵢ₋₁)*
  end
  append!(eps, w[1:n - ks[end]])              # Q*
  [eps[i] - eps[i + 1] for i in 1:n - 1]
end

# =============================================================================
#  Ambient variety construction
# =============================================================================

"""
Build a PartialFlagVariety from the `ambient` field of data.json.

Each factor `[k₁,...,kₘ, n]` becomes `Fl(k₁,...,kₘ; n) = A_{n-1}` with
marked nodes `(k₁,...,kₘ)`.  Products use `ProductDynkinType`.
"""
function build_ambient(factors::Vector)
  # Parse each factor
  factor_types = []
  factor_marks = Vector{Int}[]
  for f in factors
    fv = Int.(f)
    n = fv[end]
    ks = fv[1:end - 1]
    push!(factor_types, TypeA{n - 1})
    push!(factor_marks, ks)
  end

  if length(factor_types) == 1
    DT = factor_types[1]
    marking = Tuple(factor_marks[1])
    return partial_flag_variety(DT, marking)
  end

  # Build ProductDynkinType by nesting
  DT = factor_types[1]
  offset = rank(DT)
  all_marks = Int[factor_marks[1]...]

  for i in 2:length(factor_types)
    DT = ProductDynkinType{Tuple{DT,factor_types[i]}}
    append!(all_marks, [m + offset for m in factor_marks[i]])
    offset += rank(factor_types[i])
  end

  return partial_flag_variety(DT, Tuple(all_marks))
end

# =============================================================================
#  Bundle construction
# =============================================================================

"""
Build a `CompletelyReducibleBundle` from the `bundle` field and `ambient` factors.

The bundle field encodes:
  bundle[1]           = the bundle (always length 1 wrapper)
  bundle[1][i]        = i-th direct summand (a box product across factors)
  bundle[1][i][j]     = weight contribution from factor j (length-1 list)
  bundle[1][i][j][1]  = GL(nⱼ) weight vector
"""
function build_bundle(X, factors::Vector, bundle_data::Vector)
  MDT = marked_type(X)
  DT = PartialFlagVarieties._ambient_type(MDT)

  # Parse factor dimensions and step sizes
  factor_ks = Vector{Int}[]
  factor_ns = Int[]
  factor_ranks = Int[]
  for f in factors
    fv = Int.(f)
    n = fv[end]
    ks = fv[1:end - 1]
    push!(factor_ks, ks)
    push!(factor_ns, n)
    push!(factor_ranks, n - 1)
  end

  summands = IrrepLevi{MDT}[]

  # bundle_data[1] is the direct sum
  isempty(bundle_data) && return CompletelyReducibleBundle{MDT}(X, summands)
  ds = bundle_data[1]
  (ds isa Vector && isempty(ds)) && return CompletelyReducibleBundle{MDT}(X, summands)

  for summand in ds
    # summand[j] = weight contribution from factor j
    # Each summand[j] is a length-1 list containing the weight vector
    omega_coords = Int[]
    for (j, factor_weight) in enumerate(summand)
      w = Int.(factor_weight[1])
      ω = gl_weight_to_omega_flag(factor_ks[j], factor_ns[j], w)
      append!(omega_coords, ω)
    end

    λ = WeightLatticeElem(DT, omega_coords)
    push!(summands, IrrepLevi(MDT, λ))
  end

  CompletelyReducibleBundle{MDT}(X, summands)
end

# =============================================================================
#  Hodge diamond comparison
# =============================================================================

"""
Parse the reference Hodge diamond from data.json format.

Returns a 5×5 Matrix{BigInt} indexed as `[p+1, q+1]` = h^{p,q}.
Rows 0–2 have 5 entries; rows 3–4 have only 3, filled by Serre duality.
"""
function parse_reference_hodge(hodge_data::Vector)
  H = zeros(BigInt, 5, 5)
  known = falses(5, 5)
  for (q, row) in enumerate(hodge_data)
    for (p, val) in enumerate(row)
      if val isa Number
        H[p, q] = BigInt(round(Int, val))
        known[p, q] = true
      end
    end
  end
  # Fill unknown entries via Serre duality: h^{p,q} = h^{d-p,d-q}
  for p in 0:4
    for q in 0:4
      pi, qi = p + 1, q + 1
      sp, sq = 4 - p + 1, 4 - q + 1
      if !known[pi, qi] && sp >= 1 && sq >= 1 && known[sp, sq]
        H[pi, qi] = H[sp, sq]
        known[pi, qi] = true
      end
    end
  end
  H
end

"""Compare computed and reference Hodge diamonds; return number of mismatches."""
function compare_hodge(computed::Matrix{BigInt}, reference::Matrix{BigInt})
  mismatches = 0
  for p in 0:4
    for q in 0:4
      if computed[p + 1, q + 1] != reference[p + 1, q + 1]
        mismatches += 1
      end
    end
  end
  mismatches
end

"""Check if the reference Hodge data contains symbolic (non-numeric) values."""
function has_symbolic_hodge(hodge_data::Vector)
  for row in hodge_data
    for val in row
      val isa AbstractString && return true
    end
  end
  false
end

# =============================================================================
#  Ambient Hodge numbers (for empty bundle case)
# =============================================================================

"""
Hodge numbers of a flag variety G/P (or product thereof).

For rational homogeneous varieties:
  h^{p,q} = 0 for p ≠ q
  h^{p,p} = b_{2p}  (Betti numbers are concentrated in even degrees)

`betti_numbers(X)` returns the (d+1)-element vector [b₀, b₂, b₄, ...].
"""
function ambient_hodge_numbers(X)
  d = dimension(X)
  d == 4 || error("ambient dimension $d ≠ 4")
  betti = betti_numbers(X)
  H = zeros(BigInt, 5, 5)
  for p in 0:4
    H[p + 1, p + 1] = betti[p + 1]  # betti[p+1] = b_{2p}
  end
  H
end

# =============================================================================
#  Main computation
# =============================================================================

function process_entry(entry, idx)
  factors = entry["ambient"]
  bundle_data = entry["bundle"]

  # Skip entries with symbolic Hodge values
  if has_symbolic_hodge(entry["hodge"])
    return (status=:symbolic, idx=idx, msg="symbolic Hodge values")
  end

  ref_hodge = parse_reference_hodge(entry["hodge"])

  # Build ambient
  X = build_ambient(factors)
  d = dimension(X)

  # Empty bundle: ambient is the 4-fold
  is_empty_bundle = (bundle_data isa Vector && length(bundle_data) == 1 &&
                     bundle_data[1] isa Vector && isempty(bundle_data[1]))

  if is_empty_bundle
    d == 4 || return (status=:dim_error, idx=idx, msg="ambient dim $d ≠ 4")
    H = ambient_hodge_numbers(X)
    mismatches = compare_hodge(H, ref_hodge)
    return (status=:ok, idx=idx, hodge=H, ref=ref_hodge, mismatches=mismatches,
      h11=H[2, 2], h22=H[3, 3], h13=H[2, 4])
  end

  # Build bundle
  E = build_bundle(X, factors, bundle_data)
  r = Int(rank_bundle(E))

  if d - r != 4
    return (status=:dim_error, idx=idx, msg="dim(Z) = $(d-r) ≠ 4")
  end

  Z = zero_locus(E)

  # Compute Hodge numbers
  H = hodge_numbers(Z)
  mismatches = compare_hodge(H, ref_hodge)
  return (status=:ok, idx=idx, hodge=H, ref=ref_hodge, mismatches=mismatches,
    h11=H[2, 2], h22=H[3, 3], h13=H[2, 4])
end

# =============================================================================
#  Main entry point
# =============================================================================

function main(; datafile=nothing, max_entries=nothing, max_rank=nothing)
  # Find data.json
  if datafile === nothing
    candidates = [
      joinpath(@__DIR__, "..", "..", "fano-fourfolds", "data.json"),
      joinpath(homedir(), "Documents", "Projects", "fano-fourfolds", "data.json"),
    ]
    for c in candidates
      if isfile(c)
        datafile = c
        break
      end
    end
    datafile === nothing && error("data.json not found. Pass path as argument.")
  end

  println("\nLoading $datafile...")
  data = JSON.parsefile(datafile)
  n_total = length(data)
  println("  $n_total entries loaded.\n")

  # Sort by bundle rank (proxy for complexity), then truncate
  bundle_ranks = map(data) do entry
    bd = entry["bundle"]
    (bd isa Vector && length(bd) == 1 && bd[1] isa Vector && isempty(bd[1])) && return 0
    isempty(bd) && return 0
    ds = bd[1]
    return length(ds)
  end
  order = sortperm(bundle_ranks)
  data = data[order]

  # Filter by max bundle rank
  if max_rank !== nothing
    data = filter(data) do entry
      bd = entry["bundle"]
      (bd isa Vector && length(bd) == 1 && bd[1] isa Vector && isempty(bd[1])) && return true
      isempty(bd) && return true
      return length(bd[1]) <= max_rank
    end
  end

  # Optional entry count limit (applied AFTER sorting)
  if max_entries !== nothing
    data = data[1:min(max_entries, length(data))]
  end

  n_run = length(data)
  println("Processing $n_run entries (sorted by bundle rank)...\n")

  results_ok = []
  results_fail = []
  results_mismatch = []
  results_symbolic = []
  n_processed = 0

  p = Progress(n_run; dt=1.0, showspeed=true, enabled=stderr isa Base.TTY)

  for (i, entry) in enumerate(data)
    try
      res = process_entry(entry, i)
      if res.status == :ok
        n_processed += 1
        if res.mismatches > 0
          push!(results_mismatch, res)
        else
          push!(results_ok, res)
        end
      elseif res.status == :symbolic
        push!(results_symbolic, res)
      else
        push!(results_fail, (status=:skip, idx=i, msg=res.msg))
      end
    catch e
      msg = sprint(showerror, e)
      push!(results_fail, (status=:error, idx=i, msg=first(msg, 120)))
    end

    next!(p; showvalues=[
      (:processed, n_processed),
      (:ok, length(results_ok)),
      (:mismatch, length(results_mismatch)),
      (:symbolic, length(results_symbolic)),
      (:failed, length(results_fail)),
    ])
  end

  finish!(p)

  # Summary
  println("\n" * "="^60)
  println("  RESULTS")
  println("="^60)
  @printf("  Total entries: %d\n", n_run)
  @printf("  Processed:     %d\n", n_processed)
  @printf("  Matched:       %d  ✓\n", length(results_ok))
  @printf("  Mismatched:    %d  ✗\n", length(results_mismatch))
  @printf("  Symbolic:      %d  (skipped — parametric Hodge numbers)\n", length(results_symbolic))
  @printf("  Skipped/Error: %d\n", length(results_fail))

  if !isempty(results_mismatch)
    println("\nMismatched entries:")
    for r in results_mismatch[1:min(10, end)]
      println("  Entry $(r.idx): $(r.mismatches) Hodge number mismatches")
      println("    computed h¹¹=$(r.h11) h²²=$(r.h22) h¹³=$(r.h13)")
      println("    expected h¹¹=$(r.ref[2,2]) h²²=$(r.ref[3,3]) h¹³=$(r.ref[2,4])")
    end
  end

  if !isempty(results_fail)
    println("\nFailed entries (first 10):")
    for r in results_fail[1:min(10, end)]
      println("  Entry $(r.idx) [$(r.status)]: $(r.msg)")
    end
  end

  println()
end

# Run with defaults or parse CLI args
if abspath(PROGRAM_FILE) == @__FILE__
  datafile = length(ARGS) >= 1 ? ARGS[1] : nothing
  max_entries = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing
  max_rank = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : nothing
  main(; datafile, max_entries, max_rank)
end
