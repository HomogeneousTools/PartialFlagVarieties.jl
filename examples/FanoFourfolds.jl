# ═══════════════════════════════════════════════════════════════════════════════
#  FanoFourfolds.jl — Hodge numbers of Fano fourfolds from data.json
#
#  Reads the dataset of Fano fourfold zero loci in products of type A
#  flag varieties (from Fatighenti–Mongardi et al.), constructs each
#  ambient variety and equivariant bundle in PartialFlagVarieties.jl,
#  computes Hodge numbers via the Koszul complex, and compares with the
#  reference values.
#
#  For entries with indeterminate Hodge numbers, the computation returns
#  symbolic AffineExpr values (x_0, x_1, ...) which are compared with the
#  symbolic expressions in the reference data.
#
#  Uses Distributed for parallel processing with pmap.
#
#  Usage:
#    julia --project=. examples/FanoFourfolds.jl [data.json path] [max_entries]
# ═══════════════════════════════════════════════════════════════════════════════

using Distributed

# Add worker processes before loading packages on them
const NWORKERS = parse(Int, get(ENV, "JULIA_NUM_PROCS", "8"))
if nworkers() < NWORKERS
  addprocs(NWORKERS - nworkers() + 1)
end

@everywhere begin

using PartialFlagVarieties
using Lie
using JSON

# =============================================================================
#  GL(n) weight -> omega coordinates (general flag)
# =============================================================================

"""
Convert a Schubert2-convention GL(n) weight vector to omega coordinates.
"""
function gl_weight_to_omega_flag(ks::Vector{Int}, n::Int, w::Vector{Int})
  eps = Int[]
  append!(eps, w[n - ks[1] + 1:n])
  for i in 2:length(ks)
    append!(eps, w[n - ks[i] + 1:n - ks[i - 1]])
  end
  append!(eps, w[1:n - ks[end]])
  [eps[i] - eps[i + 1] for i in 1:n - 1]
end

# =============================================================================
#  Ambient variety construction
# =============================================================================

"""
Build a PartialFlagVariety from the `ambient` field of data.json.
"""
function build_ambient(factors::Vector)
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
"""
function build_bundle(X, factors::Vector, bundle_data::Vector)
  mdt = marked_dynkin_type(X)
  DT = PartialFlagVarieties._ambient_type(mdt)

  factor_ks = Vector{Int}[]
  factor_ns = Int[]
  for f in factors
    fv = Int.(f)
    push!(factor_ks, fv[1:end - 1])
    push!(factor_ns, fv[end])
  end

  summands = IrrepLevi[]
  isempty(bundle_data) && return CompletelyReducibleBundle(X, summands)
  ds = bundle_data[1]
  (ds isa Vector && isempty(ds)) && return CompletelyReducibleBundle(X, summands)

  for summand in ds
    omega_coords = Int[]
    for (j, factor_weight) in enumerate(summand)
      w = Int.(factor_weight[1])
      omega = gl_weight_to_omega_flag(factor_ks[j], factor_ns[j], w)
      append!(omega_coords, omega)
    end
    lam = WeightLatticeElem(DT, omega_coords)
    push!(summands, IrrepLevi(mdt, lam))
  end

  CompletelyReducibleBundle(X, summands)
end

# =============================================================================
#  Parsing symbolic reference values
# =============================================================================

"""
Parse a symbolic expression string from data.json (e.g. "2 * x_0 + 16")
into an AffineExpr.
"""
function parse_symbolic_expr(s::AbstractString)
  s = strip(s)
  constant = BigInt(0)
  coeffs = Dict{Int,BigInt}()

  # Tokenize: split on + or - (keeping the sign)
  tokens = String[]
  current = ""
  for (i, c) in enumerate(s)
    if c in ('+', '-') && i > 1 && !isempty(strip(current))
      push!(tokens, strip(current))
      current = string(c)
    else
      current *= c
    end
  end
  !isempty(strip(current)) && push!(tokens, strip(current))

  for tok in tokens
    tok = strip(tok)
    isempty(tok) && continue

    # Strip leading '+' (kept by tokenizer for positive terms after split)
    if startswith(tok, "+")
      tok = strip(tok[2:end])
      isempty(tok) && continue
    end

    # "2 * x_0", "-2 * x_0"
    m = match(r"^(-?\s*\d+)\s*\*\s*x_(\d+)$", tok)
    if m !== nothing
      coeff = parse(BigInt, replace(m.captures[1], " " => ""))
      var_id = parse(Int, m.captures[2])
      coeffs[var_id] = get(coeffs, var_id, BigInt(0)) + coeff
      continue
    end

    # "x_0", "-x_0"
    m = match(r"^(-?\s*)x_(\d+)$", tok)
    if m !== nothing
      sign_str = strip(m.captures[1])
      coeff = sign_str == "-" ? BigInt(-1) : BigInt(1)
      var_id = parse(Int, m.captures[2])
      coeffs[var_id] = get(coeffs, var_id, BigInt(0)) + coeff
      continue
    end

    # Plain integer
    m = match(r"^-?\s*\d+$", tok)
    if m !== nothing
      constant += parse(BigInt, replace(tok, " " => ""))
      continue
    end

    error("Cannot parse token '$tok' in expression '$s'")
  end

  filter!(p -> p.second != 0, coeffs)
  AffineExpr(constant, coeffs)
end

"""
Parse a reference hodge value: number or symbolic string -> AffineExpr.
"""
function parse_hodge_value(val)
  if val isa Number
    AffineExpr(BigInt(round(Int, val)))
  elseif val isa AbstractString
    parse_symbolic_expr(val)
  else
    error("Unexpected hodge value type: $(typeof(val))")
  end
end

"""
Parse the reference Hodge diamond from data.json format into AffineExpr matrix.
Fills by Serre duality h^{p,q} = h^{d-p,d-q} for d=4.
"""
function parse_reference_hodge(hodge_data::Vector)
  H = Matrix{AffineExpr}(undef, 5, 5)
  for i in eachindex(H)
    H[i] = AffineExpr(0)
  end
  known = falses(5, 5)

  for (q, row) in enumerate(hodge_data)
    for (p, val) in enumerate(row)
      H[p, q] = parse_hodge_value(val)
      known[p, q] = true
    end
  end

  # Fill by Serre duality: h^{p,q} = h^{d-p,d-q}
  for p in 0:4, q in 0:4
    pi, qi = p + 1, q + 1
    sp, sq = 4 - p + 1, 4 - q + 1
    if !known[pi, qi] && sp >= 1 && sq >= 1 && known[sp, sq]
      H[pi, qi] = H[sp, sq]
      known[pi, qi] = true
    end
  end

  H
end

"""Check if the reference Hodge data contains symbolic values."""
function has_symbolic_hodge(hodge_data::Vector)
  for row in hodge_data
    for val in row
      val isa AbstractString && return true
    end
  end
  false
end

# =============================================================================
#  Symbolic comparison with variable renaming
# =============================================================================

"""Rename variables in an AffineExpr according to a mapping.
Variables in `zero_vars` are set to 0."""
function rename_vars(e::AffineExpr, mapping::Dict{Int,Int};
  zero_vars::Set{Int}=Set{Int}())
  new_constant = e.constant
  new_coeffs = Dict{Int,BigInt}()
  for (var_id, coeff) in e.coeffs
    if var_id in zero_vars
      continue  # Variable specialized to 0
    elseif haskey(mapping, var_id)
      new_id = mapping[var_id]
      new_coeffs[new_id] = get(new_coeffs, new_id, BigInt(0)) + coeff
      new_coeffs[new_id] == 0 && delete!(new_coeffs, new_id)
    else
      # Unmapped variable — keep as-is
      new_coeffs[var_id] = get(new_coeffs, var_id, BigInt(0)) + coeff
      new_coeffs[var_id] == 0 && delete!(new_coeffs, var_id)
    end
  end
  AffineExpr(new_constant, new_coeffs)
end

"""
Try to find a consistent variable renaming (and zero-specialization)
from computed to reference expressions.

Each computed variable is either:
  - mapped to a reference variable (injective), or
  - specialized to 0 (the connecting map rank was actually zero).

Returns `(matches, mapping)`.
"""
function find_variable_mapping(
  computed::Matrix{AffineExpr}, reference::Matrix{AffineExpr},
)
  c_var_set = Set{Int}()
  r_var_set = Set{Int}()
  for e in computed
    for k in keys(e.coeffs)
      push!(c_var_set, k)
    end
  end
  for e in reference
    for k in keys(e.coeffs)
      push!(r_var_set, k)
    end
  end

  c_vars = sort(collect(c_var_set))
  r_vars = sort(collect(r_var_set))

  # No computed variables: direct comparison
  if isempty(c_vars)
    ok = all(computed[i] == reference[i] for i in eachindex(computed))
    return (ok, Dict{Int,Int}())
  end

  # Brute-force search over all assignments
  best_mapping = Dict{Int,Int}()
  best_zeros = Set{Int}()

  function verify(mapping, zero_set)
    for i in eachindex(computed)
      renamed = rename_vars(computed[i], mapping; zero_vars=zero_set)
      renamed != reference[i] && return false
    end
    true
  end

  function search(idx, mapping, zero_set, used_r)
    if idx > length(c_vars)
      if verify(mapping, zero_set)
        merge!(best_mapping, mapping)
        union!(best_zeros, zero_set)
        return true
      end
      return false
    end

    cv = c_vars[idx]

    # Try mapping to each unused reference variable
    for rv in r_vars
      rv in used_r && continue
      mapping[cv] = rv
      push!(used_r, rv)
      if search(idx + 1, mapping, zero_set, used_r)
        return true
      end
      delete!(mapping, cv)
      delete!(used_r, rv)
    end

    # Try setting to 0
    push!(zero_set, cv)
    if search(idx + 1, mapping, zero_set, used_r)
      return true
    end
    delete!(zero_set, cv)

    false
  end

  found = search(1, Dict{Int,Int}(), Set{Int}(), Set{Int}())
  (found, best_mapping)
end

# =============================================================================
#  Ambient Hodge numbers (for empty bundle case)
# =============================================================================

function ambient_hodge_numbers(X)
  d = dimension(X)
  d == 4 || error("ambient dimension $d != 4")
  betti = betti_numbers(X)
  H = Matrix{AffineExpr}(undef, 5, 5)
  for i in eachindex(H)
    H[i] = AffineExpr(0)
  end
  for p in 0:4
    H[p + 1, p + 1] = AffineExpr(betti[p + 1])
  end
  H
end

# =============================================================================
#  Main computation
# =============================================================================

function process_entry(entry_with_idx)
  entry, idx = entry_with_idx
  factors = entry["ambient"]
  bundle_data = entry["bundle"]
  is_sym = has_symbolic_hodge(entry["hodge"])
  ref_hodge = parse_reference_hodge(entry["hodge"])

  X = build_ambient(factors)
  d = dimension(X)

  is_empty_bundle = (bundle_data isa Vector && length(bundle_data) == 1 &&
                     bundle_data[1] isa Vector && isempty(bundle_data[1]))

  if is_empty_bundle
    d == 4 || return (status=:dim_error, idx=idx, msg="ambient dim $d != 4",
      is_symbolic=is_sym, hodge=nothing, ref=ref_hodge,
      matches=false, mapping=Dict{Int,Int}())
    H = ambient_hodge_numbers(X)
    (matches, mapping) = find_variable_mapping(H, ref_hodge)
    return (status=:ok, idx=idx, hodge=H, ref=ref_hodge,
      matches=matches, mapping=mapping, is_symbolic=is_sym)
  end

  E = build_bundle(X, factors, bundle_data)
  r = Int(rank_bundle(E))

  if d - r != 4
    return (status=:dim_error, idx=idx, msg="dim(Z) = $(d-r) != 4",
      is_symbolic=is_sym, hodge=nothing, ref=ref_hodge,
      matches=false, mapping=Dict{Int,Int}())
  end

  Z = zero_locus(E)
  H = hodge_numbers_symbolic(Z)
  (matches, mapping) = find_variable_mapping(H, ref_hodge)

  return (status=:ok, idx=idx, hodge=H, ref=ref_hodge,
    matches=matches, mapping=mapping, is_symbolic=is_sym)
end

end  # @everywhere

# =============================================================================
#  Pretty-print helpers
# =============================================================================

using Printf

function format_hodge_matrix(H::Matrix{AffineExpr})
  d = size(H, 1) - 1
  lines = String[]
  for q in 0:d
    entries = String[]
    for p in 0:d
      push!(entries, sprint(show, H[p + 1, q + 1]))
    end
    push!(lines, join(entries, "  "))
  end
  join(lines, "\n")
end

function print_discrepancy(r; io=stdout)
  println(io, "  Entry $(r.idx): $(r.is_symbolic ? "symbolic" : "numeric")")
  if r.hodge !== nothing
    println(io, "  Computed:")
    println(io, "    ", replace(format_hodge_matrix(r.hodge), "\n" => "\n    "))
  end
  println(io, "  Reference:")
  println(io, "    ", replace(format_hodge_matrix(r.ref), "\n" => "\n    "))
  println(io)
end

# =============================================================================
#  Main entry point
# =============================================================================

function main(; datafile=nothing, max_entries=500)
  if datafile === nothing
    candidates = [
      joinpath(@__DIR__, "FanoFourfolds.json"),
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
  println("  $n_total entries loaded.")

  # Sort by bundle rank (proxy for codimension / complexity)
  bundle_ranks = map(data) do entry
    bd = entry["bundle"]
    (bd isa Vector && length(bd) == 1 && bd[1] isa Vector && isempty(bd[1])) && return 0
    isempty(bd) && return 0
    return length(bd[1])
  end
  order = sortperm(bundle_ranks)
  data = data[order]

  # Trim to max_entries
  data = data[1:min(max_entries, end)]
  n_run = length(data)
  println("  Processing $n_run entries (sorted by bundle rank).")
  println("  Using $(nworkers()) worker process(es).\n")

  # Process in parallel with pmap
  entries_with_indices = collect(zip(data, 1:n_run))
  results = pmap(process_entry, entries_with_indices)

  # Categorize results
  ok_numeric = []
  ok_symbolic = []
  mismatch_numeric = []
  mismatch_symbolic = []
  failed = []

  for r in results
    if r.status == :ok
      if r.matches
        r.is_symbolic ? push!(ok_symbolic, r) : push!(ok_numeric, r)
      else
        r.is_symbolic ? push!(mismatch_symbolic, r) : push!(mismatch_numeric, r)
      end
    else
      push!(failed, r)
    end
  end

  # Summary
  println("=" ^ 60)
  println("  RESULTS")
  println("=" ^ 60)
  @printf("  Total entries processed: %d\n", n_run)
  @printf("  Numeric matches:    %4d\n", length(ok_numeric))
  @printf("  Symbolic matches:   %4d  (with variable renaming)\n", length(ok_symbolic))
  @printf("  Numeric mismatches: %4d\n", length(mismatch_numeric))
  @printf("  Symbolic mismatches:%4d\n", length(mismatch_symbolic))
  @printf("  Errors/skipped:     %4d\n", length(failed))
  println("=" ^ 60)

  if !isempty(mismatch_numeric)
    println("\n-- Numeric mismatches --")
    for r in mismatch_numeric
      print_discrepancy(r)
    end
  end

  if !isempty(mismatch_symbolic)
    println("\n-- Symbolic mismatches --")
    for r in mismatch_symbolic
      print_discrepancy(r)
    end
  end

  if !isempty(failed)
    println("\n-- Errors --")
    for r in failed[1:min(20, end)]
      println("  Entry $(r.idx): $(r.msg)")
    end
  end

  println()
end

# Run with defaults or parse CLI args
if abspath(PROGRAM_FILE) == @__FILE__
  datafile = length(ARGS) >= 1 ? ARGS[1] : nothing
  max_entries = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 500
  main(; datafile, max_entries)
end
