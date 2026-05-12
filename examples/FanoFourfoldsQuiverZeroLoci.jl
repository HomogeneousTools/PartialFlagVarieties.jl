# ═══════════════════════════════════════════════════════════════════════════════
#  QuiverZeroLoci.jl — Independent recomputation of Fano fourfold quiver zero loci
#
#  Reference data from the Zenodo repository:
#    https://zenodo.org/doi/10.5281/zenodo.17742466
#
#  This script loads all 170 entries of the precompiled quiver zero loci
#  Hodge-number files (both the "raw" and the "final" versions), translates
#  them into PartialFlagVarieties.jl objects, recomputes their Hodge numbers
#  and invariants (h⁰(-K), (-K)⁴, χ(T)), and compares against the published
#  reference values.  It also computes h⁰(T) and h¹(T) where possible.
#
#  Underdetermined entries (whose Hodge numbers involve free symbolic
#  variables like z_3, z_4) are identified and flagged.
#
#  The script also compares the raw and final reference files against each
#  other to flag any inconsistencies.
#
#  Usage:
#    julia --project=. examples/QuiverZeroLoci.jl [--ids F3.10,F2.42]
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PartialFlagVarieties: AffineExpr, symbolic_variable, is_determined
using Printf
using Downloads
using PrettyTables
using ProgressMeter

# ═══════════════════════════════════════════════════════════════════════════════
#  Data structures
# ═══════════════════════════════════════════════════════════════════════════════

"""Reference data for one quiver zero locus entry."""
struct QuiverEntry
  id::String
  ambient::Vector{Vector{Int}}
  bundle::Vector{Vector{Vector{Int}}}
  hodge::Matrix{AffineExpr}
  invariants::Vector{Int}  # [h⁰(-K), (-K)⁴, χ(T)]
end

"""Result of recomputing one entry."""
struct ComputationResult
  id::String
  ambient_desc::String
  hodge_computed::Union{Nothing,Matrix{AffineExpr}}
  hodge_reference::Matrix{AffineExpr}
  hodge_match::Bool
  underdetermined::Bool
  h0_antiK::Union{Nothing,Int}
  antiK4::Union{Nothing,Int}
  h0_T::Union{Nothing,Int}
  h1_T::Union{Nothing,Int}
  chi_T::Union{Nothing,Int}
  # Hodge numbers h^{p,q} as display strings ("?" if not computed)
  h11::String
  h12::String
  h13::String
  h22::String
  ref_h0_antiK::Int
  ref_antiK4::Int
  ref_chi_T::Int
  inv_match::Bool
  elapsed::Float64
  error_msg::Union{Nothing,String}
end

# ═══════════════════════════════════════════════════════════════════════════════
#  M2 parsing utilities
# ═══════════════════════════════════════════════════════════════════════════════

"""Skip whitespace in `s` starting from byte index `i`."""
function _skip_ws(s::AbstractString, i::Int)
  while i <= lastindex(s) && isspace(s[i])
    i = nextind(s, i)
  end
  i
end

"""
Extract a balanced brace-delimited literal `{...}` from `s` starting at
position `start` (which may point at whitespace before the opening `{`).
Returns `(substring, next_index)`.
"""
function _extract_braced(s::AbstractString, start::Int)
  i = _skip_ws(s, start)
  i <= lastindex(s) && s[i] == '{' ||
    error("Expected '{' near position $i: …$(s[max(1, i - 10):min(end, i + 20)])…")
  depth = 0
  j = i
  while j <= lastindex(s)
    ch = s[j]
    if ch == '{'
      depth += 1
    elseif ch == '}'
      depth -= 1
      depth == 0 && return s[i:j], nextind(s, j)
    end
    j = nextind(s, j)
  end
  error("Unterminated brace starting at position $i")
end

"""
Parse a nested M2 list literal like `{{1, 2}, {3, 4}}` into nested vectors,
using `atom_parser` for leaf values.  Missing values (empty fields between
commas, as in `{0,1,,}`) are represented as `missing`.
"""
function _parse_m2_list(s::AbstractString, atom_parser::Function)
  s = strip(s)
  isempty(s) && error("Empty M2 list literal")
  s[1] == '{' && s[end] == '}' || error("Expected braced list: $(first(s, 60))")
  inner = s[nextind(s, 1):prevind(s, lastindex(s))]

  items = Any[]
  depth = 0
  token_start = firstindex(inner)
  i = firstindex(inner)
  last_idx = lastindex(inner)

  while i <= last_idx
    ch = inner[i]
    if ch == '{'
      depth += 1
    elseif ch == '}'
      depth -= 1
    elseif ch == ',' && depth == 0
      token = strip(inner[token_start:prevind(inner, i)])
      if isempty(token)
        push!(items, missing)
      elseif startswith(token, "{")
        push!(items, _parse_m2_list(token, atom_parser))
      else
        push!(items, atom_parser(token))
      end
      token_start = nextind(inner, i)
      i = nextind(inner, i)
      continue
    end
    i = nextind(inner, i)
  end

  # Handle last field
  if token_start <= last_idx
    token = strip(inner[token_start:last_idx])
    if isempty(token)
      push!(items, missing)
    elseif startswith(token, "{")
      push!(items, _parse_m2_list(token, atom_parser))
    else
      push!(items, atom_parser(token))
    end
  elseif !isempty(inner) && strip(inner[end:end]) == ","
    # Trailing comma with nothing after → missing
    push!(items, missing)
  elseif !isempty(items)
    # Trailing comma already processed → push missing
    # (this catches the case where the inner string ends right after a comma)
    last_char_idx = lastindex(inner)
    if last_char_idx >= 1 && inner[last_char_idx] == ','
      push!(items, missing)
    end
  end

  items
end

_parse_int(token::AbstractString) = parse(Int, strip(token))

"""
Parse an affine expression token like `z_3+3`, `2*z_3+8`, or a plain integer.
Returns an `AffineExpr`.  Empty or whitespace-only tokens yield `missing`.
"""
function _parse_affine(token::AbstractString)
  stripped = replace(strip(token), " " => "")
  isempty(stripped) && return missing

  v = tryparse(BigInt, stripped)
  v !== nothing && return AffineExpr(v)

  constant = BigInt(0)
  coeffs = Dict{Int,BigInt}()
  i = firstindex(stripped)
  last_idx = lastindex(stripped)

  while i <= last_idx
    sign = BigInt(1)
    if stripped[i] == '+'
      i = nextind(stripped, i)
    elseif stripped[i] == '-'
      sign = BigInt(-1)
      i = nextind(stripped, i)
    end

    start = i
    while i <= last_idx && stripped[i] != '+' && stripped[i] != '-'
      i = nextind(stripped, i)
    end
    term = stripped[start:prevind(stripped, i)]
    isempty(term) && continue

    var_match = match(r"^(?:(\d+)\*)?([A-Za-z]+)_(\d+)$", term)
    if var_match !== nothing
      coeff =
        var_match.captures[1] === nothing ? BigInt(1) : parse(BigInt, var_match.captures[1])
      var_id = parse(Int, var_match.captures[3])
      coeffs[var_id] = get(coeffs, var_id, BigInt(0)) + sign * coeff
      coeffs[var_id] == 0 && delete!(coeffs, var_id)
      continue
    end

    if occursin(r"^-?\d+$", term)
      constant += sign * parse(BigInt, term)
      continue
    end

    error("Cannot parse affine token: '$token' (term='$term')")
  end

  AffineExpr(constant, coeffs)
end

"""
Complete a partial Hodge diamond by Hodge symmetry h^{p,q} = h^{q,p}
and Serre duality h^{p,q} = h^{d-p,d-q}.
"""
function _complete_hodge_matrix(raw_rows)
  d = length(raw_rows) - 1
  H = Matrix{Union{Missing,AffineExpr}}(missing, d + 1, d + 1)

  for p in 0:d
    row = raw_rows[p + 1]
    for q in 0:min(d, length(row) - 1)
      v = row[q + 1]
      if v !== missing
        H[p + 1, q + 1] = v isa Integer ? AffineExpr(v) : v
      end
    end
  end

  changed = true
  while changed
    changed = false
    for p in 0:d, q in 0:d
      for (pp, qq) in ((q, p), (d - p, d - q), (d - q, d - p))
        0 <= pp <= d && 0 <= qq <= d || continue
        a = H[p + 1, q + 1]
        b = H[pp + 1, qq + 1]
        if a === missing && b !== missing
          H[p + 1, q + 1] = b
          changed = true
        end
      end
    end
  end

  any(ismissing, H) && error("Could not complete Hodge matrix from symmetry data")
  Matrix{AffineExpr}([H[i, j]::AffineExpr for i in 1:(d + 1), j in 1:(d + 1)])
end

# ═══════════════════════════════════════════════════════════════════════════════
#  File parsers
# ═══════════════════════════════════════════════════════════════════════════════

"""
Find the braced value of a variable assignment `VAR = {...}` in `text`.
Returns `(braced_string, next_position)`.
"""
function _find_assignment(text::AbstractString, var::AbstractString)
  pattern = Regex("\\b$(var)\\s*=\\s*")
  m = match(pattern, text)
  m === nothing && error("Assignment '$var = {...}' not found")
  _extract_braced(text, m.offset + lastindex(m.match))
end

"""
Parse a text block's variable assignments (X, F, HN, invs) into a `QuiverEntry`.
"""
function _parse_block(id::AbstractString, block::AbstractString)
  x_str, _ = _find_assignment(block, "X")
  f_str, _ = _find_assignment(block, "F")
  hn_str, _ = _find_assignment(block, "HN")
  invs_str, _ = _find_assignment(block, "invs")

  ambient_raw = _parse_m2_list(x_str, _parse_int)
  ambient = if !isempty(ambient_raw) && ambient_raw[1] isa AbstractVector
    [Int.(v) for v in ambient_raw]
  else
    [Int.(ambient_raw)]
  end

  bundle_raw = _parse_m2_list(f_str, _parse_int)
  bundle = Vector{Vector{Int}}[]
  for summand in bundle_raw
    if summand isa AbstractVector && !isempty(summand) && summand[1] isa AbstractVector
      push!(bundle, [Int.(w) for w in summand])
    else
      push!(bundle, [Int.(summand)])
    end
  end

  hodge_rows = _parse_m2_list(hn_str, _parse_affine)
  hodge = _complete_hodge_matrix(hodge_rows)

  invs = Int.(_parse_m2_list(invs_str, _parse_int))

  QuiverEntry(id, ambient, bundle, hodge, invs)
end

"""
    load_final_file(path) -> Dict{String, QuiverEntry}

Parse the multi-line `precompiledQuiverZeroLociHodgeNumbers` file.
Each entry consists of `X = ..., F = ..., HN = ..., invs = ...` blocks
followed by `quiverZeroLoci#"ID"={X,F,HN,invs};`.
"""
function load_final_file(text::AbstractString)
  flat = replace(text, '\n' => ' ')

  entries = Dict{String,QuiverEntry}()
  marker_re = r"quiverZeroLoci#\"([^\"]+)\"=\{X,F,HN,invs\};"
  markers = collect(eachmatch(marker_re, flat))

  prev_end = 1
  for m in markers
    id = m.captures[1]
    block = flat[prev_end:(m.offset - 1)]
    prev_end = m.offset + lastindex(m.match)
    try
      entries[id] = _parse_block(id, block)
    catch err
      @warn "Failed to parse final entry $id" exception = (err, catch_backtrace())
    end
  end

  entries
end

"""
    load_raw_file(text) -> Dict{String, QuiverEntry}

Parse the compact `precompiledQuiverZeroLociHodgeNumbersRaw` content.
Entries are top-level brace groups `{"ID", {{ambient},{bundle}}, {hodge}, {invs}}`.
"""
function load_raw_file(text::AbstractString)
  flat = replace(text, '\n' => ' ')
  entries = Dict{String,QuiverEntry}()

  i = _skip_ws(flat, 1)
  while i <= lastindex(flat)
    if flat[i] == '{'
      group, next_i = _extract_braced(flat, i)
      try
        _parse_raw_entry!(entries, group)
      catch err
        snippet = first(group, 60)
        @warn "Failed to parse raw entry starting with: $snippet" exception = (
          err, catch_backtrace()
        )
      end
      i = _skip_ws(flat, next_i)
      if i <= lastindex(flat) && flat[i] == ','
        i = _skip_ws(flat, nextind(flat, i))
      end
    else
      i = nextind(flat, i)
    end
  end

  entries
end

"""Parse one raw-format entry and insert into `entries`."""
function _parse_raw_entry!(entries::Dict{String,QuiverEntry}, group::AbstractString)
  inner = String(strip(group))
  inner = inner[nextind(inner, 1):prevind(inner, lastindex(inner))]

  m = match(r"^\s*\"([^\"]+)\"", inner)
  m === nothing && error("No quoted ID in raw entry: $(first(inner, 40))")
  id = m.captures[1]
  pos = _skip_ws(inner, m.offset + lastindex(m.match))

  inner[pos] == ',' || error("Expected comma after ID in entry $id")
  pos = _skip_ws(inner, nextind(inner, pos))

  spec_str, pos = _extract_braced(inner, pos)
  pos = _skip_ws(inner, pos)
  if pos <= lastindex(inner) && inner[pos] == ','
    pos = _skip_ws(inner, nextind(inner, pos))
  end

  hodge_str, pos = _extract_braced(inner, pos)
  pos = _skip_ws(inner, pos)
  if pos <= lastindex(inner) && inner[pos] == ','
    pos = _skip_ws(inner, nextind(inner, pos))
  end

  invs_str, _ = _extract_braced(inner, pos)

  # Parse spec = {ambient, bundle}
  spec_inner = spec_str[nextind(spec_str, 1):prevind(spec_str, lastindex(spec_str))]
  amb_str, sp = _extract_braced(spec_inner, _skip_ws(spec_inner, 1))
  sp = _skip_ws(spec_inner, sp)
  if sp <= lastindex(spec_inner) && spec_inner[sp] == ','
    sp = _skip_ws(spec_inner, nextind(spec_inner, sp))
  end
  bun_str, _ = _extract_braced(spec_inner, sp)

  ambient_raw = _parse_m2_list(amb_str, _parse_int)
  ambient = if !isempty(ambient_raw) && ambient_raw[1] isa AbstractVector
    [Int.(v) for v in ambient_raw]
  else
    [Int.(ambient_raw)]
  end

  bundle_raw = _parse_m2_list(bun_str, _parse_int)
  bundle = Vector{Vector{Int}}[]
  for summand in bundle_raw
    if summand isa AbstractVector && !isempty(summand) && summand[1] isa AbstractVector
      push!(bundle, [Int.(w) for w in summand])
    else
      push!(bundle, [Int.(summand)])
    end
  end

  hodge_rows = _parse_m2_list(hodge_str, _parse_affine)
  hodge = _complete_hodge_matrix(hodge_rows)

  invs = Int.(_parse_m2_list(invs_str, _parse_int))

  entries[id] = QuiverEntry(id, ambient, bundle, hodge, invs)
  nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Translation: GL(n) encoding → PartialFlagVarieties.jl types
# ═══════════════════════════════════════════════════════════════════════════════

"""
Convert a GL(n) weight vector (standard basis ε₁,…,εₙ) for Gr(k,n)
into omega-coordinates (fundamental weight basis) for type A_{n-1}.

The Schubert2 convention places the tautological subbundle in positions
n-k+1 through n, so we permute accordingly.
"""
function gl_to_omega(k::Int, n::Int, w::Vector{Int})
  length(w) == n || error("Expected weight of length $n, got $(length(w))")
  eps = vcat(w[(n - k + 1):n], w[1:(n - k)])
  [eps[i] - eps[i + 1] for i in 1:(n - 1)]
end

"""
Build a `PartialFlagVariety` from Grassmannian factors `[[k₁,n₁], [k₂,n₂], ...]`.
"""
function build_ambient(factors::Vector{Vector{Int}})
  length(factors) >= 1 || error("At least one Grassmannian factor required")

  factor_types = DataType[]
  factor_marks = Vector{Int}[]
  for (k, n) in factors
    1 <= k < n || error("Invalid Grassmannian Gr($k, $n)")
    push!(factor_types, TypeA{n - 1})
    push!(factor_marks, [k])
  end

  if length(factor_types) == 1
    return partial_flag_variety(factor_types[1], factor_marks[1])
  end

  DT = factor_types[1]
  offset = rank(DT)
  marked = Int[factor_marks[1]...]

  for i in 2:length(factor_types)
    DT = ProductDynkinType{Tuple{DT,factor_types[i]}}
    append!(marked, [m + offset for m in factor_marks[i]])
    offset += rank(factor_types[i])
  end

  partial_flag_variety(DT, marked)
end

"""
Build a `CompletelyReducibleBundle` on `X` from the GL-weight encoding
of bundle summands, using the coefficient-list `CompletelyReducibleBundle`
constructor directly.
"""
function build_bundle(
  X::PartialFlagVariety,
  factors::Vector{Vector{Int}},
  bundle_data::Vector{Vector{Vector{Int}}},
)
  summands = Vector{Int}[]

  for summand in bundle_data
    length(summand) == length(factors) || error(
      "Bundle summand has $(length(summand)) factors, expected $(length(factors))"
    )

    omega_coords = Int[]
    for (factor, weight) in zip(factors, summand)
      k, n = factor
      append!(omega_coords, gl_to_omega(k, n, Int.(weight)))
    end

    push!(summands, omega_coords)
  end

  isempty(summands) && error("Empty bundle data")
  CompletelyReducibleBundle(X, summands)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Invariant computation
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compute h⁰(-K_Z) = χ(-K_Z) (by Kodaira vanishing on Fano varieties).
-K_Z on X is ``(-K_X) \\otimes \\det(E)^*``.
"""
function compute_h0_antiK(Z::ZeroLocus, X::PartialFlagVariety, E::CompletelyReducibleBundle)
  antiK_Z = tensor_product(anticanonical_bundle(X), dual(det(E)))
  Int(euler_characteristic(Z, antiK_Z))
end

"""
Compute (-K_Z)⁴ via finite differences of χ(Z, (-K_Z)^{⊗t}) for t = 0,…,4.
"""
function compute_antiK_fourth(
  Z::ZeroLocus, X::PartialFlagVariety, E::CompletelyReducibleBundle
)
  L = tensor_product(anticanonical_bundle(X), dual(det(E)))

  vals = BigInt[]
  for t in 0:4
    if t == 0
      push!(vals, euler_characteristic(Z))
    else
      L_t = reduce(tensor_product, fill(L, t))
      push!(vals, euler_characteristic(Z, L_t))
    end
  end

  Int(vals[5] - 4vals[4] + 6vals[3] - 4vals[2] + vals[1])
end

"""
Compute χ(T_Z) using additivity: χ(T_Z) = χ(T_X|_Z) - χ(E|_Z).
"""
function compute_chi_tangent(
  Z::ZeroLocus, X::PartialFlagVariety, E::CompletelyReducibleBundle
)
  chi_TX = euler_characteristic(Z, tangent_bundle(X))
  chi_E = euler_characteristic(Z, E)
  Int(chi_TX - chi_E)
end

"""
Attempt to compute h⁰(T_Z) and h¹(T_Z) from the short exact sequence
  0 → T_Z → T_X|_Z → E|_Z → 0
via the associated long exact sequence in cohomology.

Returns `(h0, h1, h0_determined, h1_determined)`.
"""
function compute_tangent_cohomology(
  Z::ZeroLocus, X::PartialFlagVariety, E::CompletelyReducibleBundle
)
  coh_TX, det_TX = cohomology_on_restriction(Z, tangent_bundle(X))
  coh_E, det_E = cohomology_on_restriction(Z, E)

  if !(det_TX && det_E)
    return nothing, nothing, false, false
  end

  d = Int(dimension(Z))

  # b[i] = H^i(T_X|_Z), c[i] = H^i(E|_Z)   (i = 0, …, d)
  b = [coh_TX[i] for i in 0:d]
  c = [coh_E[i] for i in 0:d]

  # LES: 0 → a₀ → b₀ → c₀ →^{δ₀} a₁ → b₁ → c₁ → ⋯ → aₙ → bₙ → cₙ → 0
  # where a = H*(T_Z).  The connecting morphisms δᵢ satisfy
  #   aᵢ = bᵢ - cᵢ + δ_{i-1} + δᵢ
  # with δ_{-1} = δ_d = 0.  We have 0 ≤ δᵢ ≤ cᵢ.
  δ_lb = zeros(BigInt, d)
  δ_ub = [min(c[i + 1], BigInt(10_000)) for i in 0:(d - 1)]

  # Tighten bounds by propagation
  for _ in 1:20
    old_lb = copy(δ_lb)
    old_ub = copy(δ_ub)
    for i in 0:(d - 1)
      δ_prev_lb = i == 0 ? BigInt(0) : δ_lb[i]
      δ_prev_ub = i == 0 ? BigInt(0) : δ_ub[i]
      # aᵢ = bᵢ - cᵢ + δ_{i-1} + δᵢ ≥ 0
      # → δᵢ ≥ cᵢ - bᵢ - δ_{i-1}  (using upper bound of δ_{i-1})
      new_lb = max(BigInt(0), c[i + 1] - b[i + 1] - δ_prev_ub)
      δ_lb[i + 1] = max(δ_lb[i + 1], new_lb)
      # Similarly: aᵢ = bᵢ - cᵢ + δ_{i-1} + δᵢ ≥ 0
      #  → δ_{i-1} ≥ cᵢ - bᵢ - δᵢ  (using upper bound of δᵢ)
    end
    for i in (d - 1):-1:0
      δ_next_lb = i == d - 1 ? BigInt(0) : δ_lb[i + 2]
      δ_next_ub = i == d - 1 ? BigInt(0) : δ_ub[i + 2]
      # a_{i+1} = b_{i+1} - c_{i+1} + δᵢ + δ_{i+1} ≥ 0
      new_lb = max(BigInt(0), c[i + 2] - b[i + 2] - δ_next_ub)
      δ_lb[i + 1] = max(δ_lb[i + 1], new_lb)
      new_ub = min(δ_ub[i + 1], b[i + 2] - c[i + 2] + δ_next_ub)
      new_ub = max(new_ub, δ_lb[i + 1])
      δ_ub[i + 1] = min(δ_ub[i + 1], new_ub)
    end
    δ_lb == old_lb && δ_ub == old_ub && break
  end

  # If uniquely determined
  if all(δ_lb[i] == δ_ub[i] for i in 1:d)
    a = BigInt[]
    for i in 0:d
      δ_prev = i == 0 ? BigInt(0) : δ_lb[i]
      δ_curr = i == d ? BigInt(0) : δ_lb[i + 1]
      push!(a, b[i + 1] - c[i + 1] + δ_prev + δ_curr)
    end
    return Int(a[1]), Int(a[2]), true, true
  end

  # Enumerate if ranges are small enough
  ranges_small = all(δ_ub[i] - δ_lb[i] <= 30 for i in 1:d)
  if !ranges_small
    return nothing, nothing, false, false
  end

  solutions = Vector{BigInt}[]
  function enumerate_deltas(idx, deltas)
    if idx > d
      a = BigInt[]
      for i in 0:d
        δ_prev = i == 0 ? BigInt(0) : deltas[i]
        δ_curr = i == d ? BigInt(0) : deltas[i + 1]
        val = b[i + 1] - c[i + 1] + δ_prev + δ_curr
        val < 0 && return nothing
        push!(a, val)
      end
      push!(solutions, a)
      return nothing
    end
    for δ in δ_lb[idx]:δ_ub[idx]
      deltas[idx] = δ
      enumerate_deltas(idx + 1, deltas)
    end
  end

  enumerate_deltas(1, zeros(BigInt, d))

  if length(solutions) == 1
    return Int(solutions[1][1]), Int(solutions[1][2]), true, true
  elseif length(solutions) > 1
    h0_vals = Set(s[1] for s in solutions)
    h1_vals = Set(s[2] for s in solutions)
    h0 = length(h0_vals) == 1 ? Int(first(h0_vals)) : nothing
    h1 = length(h1_vals) == 1 ? Int(first(h1_vals)) : nothing
    h0_det = h0 !== nothing
    h1_det = h1 !== nothing
    return h0, h1, h0_det, h1_det
  end

  nothing, nothing, false, false
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge number comparison
# ═══════════════════════════════════════════════════════════════════════════════

"""Sort key for entry IDs: `(family_number, index)`."""
function _id_key(id::AbstractString)
  m = match(r"^F(\d+)\.(\d+)$", id)
  m === nothing && return (typemax(Int), typemax(Int))
  (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
end

"""Rename symbolic variables in an `AffineExpr` according to a mapping."""
function _rename_vars(
  expr::AffineExpr,
  mapping::Dict{Int,Int};
  zero_vars::Set{Int}=Set{Int}(),
)
  constant = expr.constant
  coeffs = Dict{Int,BigInt}()
  for (var_id, coeff) in expr.coeffs
    var_id in zero_vars && continue
    new_id = get(mapping, var_id, var_id)
    coeffs[new_id] = get(coeffs, new_id, BigInt(0)) + coeff
    coeffs[new_id] == 0 && delete!(coeffs, new_id)
  end
  AffineExpr(constant, coeffs)
end

"""
Find a variable renaming bijection such that `computed` matches `reference`.
Returns `(matches::Bool, mapping::Dict{Int,Int})`.
"""
function _match_hodge(
  computed::Matrix{AffineExpr},
  reference::Matrix{AffineExpr},
)
  size(computed) == size(reference) || return false, Dict{Int,Int}()

  computed_vars = sort!(collect(Set(vcat([collect(keys(e.coeffs)) for e in computed]...))))
  reference_vars = sort!(
    collect(Set(vcat([collect(keys(e.coeffs)) for e in reference]...)))
  )

  if isempty(computed_vars) && isempty(reference_vars)
    return all(computed[i] == reference[i] for i in eachindex(computed)), Dict{Int,Int}()
  end

  best_mapping = Dict{Int,Int}()
  best_zeros = Set{Int}()
  found = Ref(false)

  function verify(mapping, zero_set)
    for i in eachindex(computed)
      _rename_vars(computed[i], mapping; zero_vars=zero_set) == reference[i] || return false
    end
    true
  end

  function search(idx, mapping, zero_set, used_ref)
    found[] && return nothing
    if idx > length(computed_vars)
      if verify(mapping, zero_set)
        merge!(best_mapping, mapping)
        union!(best_zeros, zero_set)
        found[] = true
      end
      return nothing
    end

    cv = computed_vars[idx]
    for rv in reference_vars
      rv in used_ref && continue
      mapping[cv] = rv
      push!(used_ref, rv)
      search(idx + 1, mapping, zero_set, used_ref)
      found[] && return nothing
      delete!(mapping, cv);
      delete!(used_ref, rv)
    end

    push!(zero_set, cv)
    search(idx + 1, mapping, zero_set, used_ref)
    found[] && return nothing
    delete!(zero_set, cv)
  end

  if isempty(computed_vars)
    return verify(Dict{Int,Int}(), Set{Int}()), Dict{Int,Int}()
  end

  search(1, Dict{Int,Int}(), Set{Int}(), Set{Int}())
  found[] && return true, best_mapping

  # Fallback: affine substitution check.
  # The computed variables may not biject to reference variables but there
  # may still exist an affine map  y_k = α_{k,0} + Σ_j α_{k,j} x_j  sending
  # reference variables to computed-variable expressions.  The requirement is
  # that for every matrix entry i and every "column" j ∈ {0,…,nc}, the RHS
  # vector lies in the column span of the reference-coefficient matrix.
  _check_affine_equivalence(computed, reference) && return true, Dict{Int,Int}()

  false, Dict{Int,Int}()
end

"""
Check whether two symbolic Hodge matrices are related by an affine
substitution of symbolic variables.

Given computed `c_i = a_i + Σ_j c_{i,j} x_j` and reference
`r_i = b_i + Σ_k r_{i,k} y_k`, we ask: do there exist rationals
`α_{k,j}` such that `y_k = α_{k,0} + Σ_j α_{k,j} x_j` makes all
entries match?  This is equivalent to: for each column `j`, the vector
`[c_{i,j}]_i` lies in the column span of the matrix `R` with entries
`R[i,k] = r_{i,k}`.
"""
function _check_affine_equivalence(
  computed::Matrix{AffineExpr},
  reference::Matrix{AffineExpr},
)
  size(computed) == size(reference) || return false

  n = length(computed)
  c_vars = sort!(collect(Set(k for e in computed for k in keys(e.coeffs))))
  r_vars = sort!(collect(Set(k for e in reference for k in keys(e.coeffs))))
  nc = length(c_vars)
  nr = length(r_vars)

  (nc == 0 && nr == 0) && return all(computed[i] == reference[i] for i in 1:n)
  nr == 0 && return false  # reference is fully determined but computed is not

  # Build reference coefficient matrix R (n × nr, rational)
  R = Matrix{Rational{BigInt}}(undef, n, nr)
  for i in 1:n
    for (k, rv) in enumerate(r_vars)
      R[i, k] = Rational{BigInt}(get(reference[i].coeffs, rv, BigInt(0)))
    end
  end

  rank_R = _rational_rank(R)

  # Check each "column": constant diff, then each computed variable
  for j in 0:nc
    b = Vector{Rational{BigInt}}(undef, n)
    for i in 1:n
      if j == 0
        b[i] = Rational{BigInt}(computed[i].constant - reference[i].constant)
      else
        b[i] = Rational{BigInt}(get(computed[i].coeffs, c_vars[j], BigInt(0)))
      end
    end

    # b must lie in ColSpan(R): check rank([R b]) == rank(R)
    aug = hcat(copy(R), b)
    _rational_rank(aug) == rank_R || return false
  end

  true
end

"""Compute the rank of a rational matrix via Gaussian elimination."""
function _rational_rank(M::Matrix{Rational{BigInt}})
  A = copy(M)
  rows, cols = size(A)
  r = 0
  for col in 1:cols
    pivot = findfirst(i -> A[i, col] != 0, (r + 1):rows)
    pivot === nothing && continue
    pivot += r
    r += 1
    if pivot != r
      A[r, :], A[pivot, :] = A[pivot, :], A[r, :]
    end
    for i in (r + 1):rows
      if A[i, col] != 0
        factor = A[i, col] / A[r, col]
        for c in 1:cols
          A[i, c] -= factor * A[r, c]
        end
      end
    end
  end
  r
end

"""Check whether a Hodge matrix contains underdetermined (symbolic) entries."""
is_underdetermined(H::Matrix{AffineExpr}) = any(!is_determined(e) for e in H)

# ═══════════════════════════════════════════════════════════════════════════════
#  Raw vs final comparison
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compare the raw and final reference files entry by entry.
Returns a list of `(id, field, details)` for any inconsistencies.
"""
function compare_raw_vs_final(
  final::Dict{String,QuiverEntry},
  raw::Dict{String,QuiverEntry},
)
  issues = Tuple{String,String,String}[]
  all_ids = sort!(union(collect(keys(final)), collect(keys(raw))); by=_id_key)

  for id in all_ids
    if !haskey(final, id)
      push!(issues, (id, "missing", "present in raw but not in final"))
      continue
    end
    if !haskey(raw, id)
      push!(issues, (id, "missing", "present in final but not in raw"))
      continue
    end

    f = final[id]
    r = raw[id]

    if f.ambient != r.ambient
      push!(issues, (id, "ambient", "final=$(f.ambient) raw=$(r.ambient)"))
    end
    if f.bundle != r.bundle
      push!(
        issues,
        (id, "bundle", "final=$(length(f.bundle)) summands, raw=$(length(r.bundle))"),
      )
    end
    if f.invariants != r.invariants
      push!(issues, (id, "invariants", "final=$(f.invariants) raw=$(r.invariants)"))
    end
    match_ok, _ = _match_hodge(f.hodge, r.hodge)
    if !match_ok
      push!(issues, (id, "hodge", "Hodge matrices differ even up to variable renaming"))
    end
  end

  issues
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main computation
# ═══════════════════════════════════════════════════════════════════════════════

"""Format the ambient variety as a string for display."""
_ambient_desc(factors::Vector{Vector{Int}}) = join(
  ["Gr($(k),$(n))" for (k, n) in factors], "×"
)

"""Process a single entry: build variety, compute invariants, compare."""
function process_entry(entry::QuiverEntry)
  t0 = time_ns()
  try
    X = build_ambient(entry.ambient)
    E = build_bundle(X, entry.ambient, entry.bundle)
    dim_Z = Int(dimension(X) - rank_bundle(E))
    dim_Z == 4 || error("Expected 4-fold, got dimension $dim_Z")

    Z = zero_locus(E)

    hodge = hodge_numbers_symbolic(Z)

    h0_aK = compute_h0_antiK(Z, X, E)
    aK4 = compute_antiK_fourth(Z, X, E)
    chi_T = compute_chi_tangent(Z, X, E)

    h0_T, h1_T, _, _ = compute_tangent_cohomology(Z, X, E)

    match_hodge, _ = _match_hodge(hodge, entry.hodge)

    inv_match =
      h0_aK == entry.invariants[1] &&
      aK4 == entry.invariants[2] &&
      chi_T == entry.invariants[3]

    # Extract individual Hodge numbers as strings (hodge is 1-indexed: [p+1, q+1])
    _hstr(m, p, q) = sprint(show, m[p + 1, q + 1])
    h11 = _hstr(hodge, 1, 1)
    h12 = _hstr(hodge, 1, 2)
    h13 = _hstr(hodge, 1, 3)
    h22 = _hstr(hodge, 2, 2)

    elapsed = (time_ns() - t0) / 1.0e9

    ComputationResult(
      entry.id,
      _ambient_desc(entry.ambient),
      hodge,
      entry.hodge,
      match_hodge,
      is_underdetermined(hodge),
      h0_aK, aK4, h0_T, h1_T, chi_T,
      h11, h12, h13, h22,
      entry.invariants[1], entry.invariants[2], entry.invariants[3],
      inv_match,
      elapsed,
      nothing,
    )
  catch err
    elapsed = (time_ns() - t0) / 1.0e9
    ComputationResult(
      entry.id,
      _ambient_desc(entry.ambient),
      nothing,
      entry.hodge,
      false,
      is_underdetermined(entry.hodge),
      nothing, nothing, nothing, nothing, nothing,
      "?", "?", "?", "?",
      entry.invariants[1], entry.invariants[2], entry.invariants[3],
      false,
      elapsed,
      sprint(showerror, err),
    )
  end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main entry point
# ═══════════════════════════════════════════════════════════════════════════════

# Zenodo repository: https://zenodo.org/doi/10.5281/zenodo.17742466 (DOI redirects to latest version)
const ZENODO_PRECOMPILED_COMPACT = "https://zenodo.org/records/17742467/files/precompiledQuiverZeroLociHodgeNumbers"
const ZENODO_PRECOMPILED_RAW = "https://zenodo.org/records/17742467/files/precompiledQuiverZeroLociHodgeNumbersRaw"

"""Download a file from a URL and return its contents as a string."""
function download_file(url::AbstractString, name::AbstractString)
  tmpdir = mktempdir()
  tmpfile = joinpath(tmpdir, name)
  println("  Downloading $name from Zenodo...")
  last_err = nothing
  for (attempt, timeout) in enumerate((300.0, 600.0, 900.0))
    try
      Downloads.download(url, tmpfile; timeout=timeout)
      return read(tmpfile, String)
    catch err
      last_err = err
      attempt == 3 && break
      println("    download attempt $attempt failed: $(sprint(showerror, err))")
      println("    retrying with a longer timeout...")
      sleep(attempt)
    end
  end
  throw(last_err)
end

"""
    main(; ids=nothing)

Load the quiver zero loci data from Zenodo, recompute all Hodge numbers
and invariants, and compare with the reference values.
"""
function main(;
  ids::Union{Nothing,Vector{String}}=nothing
)
  # Download data from Zenodo
  println()
  final_content = download_file(
    ZENODO_PRECOMPILED_COMPACT, "precompiledQuiverZeroLociHodgeNumbers"
  )
  raw_content = download_file(
    ZENODO_PRECOMPILED_RAW, "precompiledQuiverZeroLociHodgeNumbersRaw"
  )
  println()

  println("=" ^ 120)
  println(
    " QuiverZeroLoci.jl — Independent Recomputation of Fano Fourfold Quiver Zero Loci"
  )
  println(" Reference: https://zenodo.org/doi/10.5281/zenodo.17742466")
  println("=" ^ 120)

  # ─── Load data ─────────────────────────────────────────────────────────
  print("Loading final file... ")
  final_entries = load_final_file(final_content)
  println("$(length(final_entries)) entries")

  print("Loading raw file...   ")
  raw_entries = load_raw_file(raw_content)
  println("$(length(raw_entries)) entries")

  # ─── Compare raw vs final ──────────────────────────────────────────────
  println("\n─── Raw vs Final Comparison ───")
  issues = compare_raw_vs_final(final_entries, raw_entries)
  if isempty(issues)
    println("  No inconsistencies between raw and final files.")
  else
    println("  $(length(issues)) inconsistencies found:")
    for (id, field, detail) in issues
      println("    $id [$field]: $detail")
    end
  end

  # ─── Recomputation ─────────────────────────────────────────────────────
  all_ids = sort!(collect(keys(final_entries)); by=_id_key)
  if ids !== nothing
    all_ids = filter(id -> id in ids, all_ids)
  end

  println("\n─── Recomputation ($(length(all_ids)) entries) ───")

  t_total = time_ns()
  entries_to_process = [final_entries[id] for id in all_ids]

  results = []
  p = Progress(length(entries_to_process);
    desc="Computing entries: ",
    barlen=50,
    showspeed=false)
  for (i, entry) in enumerate(entries_to_process)
    r = process_entry(entry)
    push!(results, r)
    next!(p; showvalues=[(:entry, "$i/$(length(entries_to_process))")])
  end

  elapsed_total = (time_ns() - t_total) / 1.0e9

  # ─── Build and display results table ────────────────────────────────
  if !isempty(results)
    table_rows = []
    for r in results
      if r.error_msg !== nothing
        row = [
          r.id, r.ambient_desc, "", "", "", "", "", "", "", "", "ERROR: $(r.error_msg)"
        ]
      else
        status = if r.hodge_match && r.inv_match
          r.underdetermined ? "OK (symbolic)" : "OK"
        else
          "MISMATCH"
        end
        h0_str = r.h0_antiK === nothing ? "?" : string(r.h0_antiK)
        aK4_str = r.antiK4 === nothing ? "?" : string(r.antiK4)
        h0T_str = r.h0_T === nothing ? "?" : string(r.h0_T)
        h1T_str = r.h1_T === nothing ? "?" : string(r.h1_T)

        row = [r.id, r.ambient_desc, h0_str, aK4_str, h0T_str, h1T_str,
          r.h11, r.h12, r.h13, r.h22, status]
      end
      push!(table_rows, row)
    end

    # Convert to matrix and display with PrettyTables.
    # display_size=(-1,-1) disables terminal-size detection and prevents
    # "N rows omitted" truncation regardless of table length.
    header = ["ID", "Ambient", "h⁰(-K)", "(-K)⁴", "h⁰(T)", "h¹(T)",
      "h¹¹", "h¹²", "h¹³", "h²²", "Status"]
    table_matrix = reduce(vcat, permutedims.(table_rows))
    pretty_table(table_matrix;
      column_labels=header,
      display_size=(-1, -1),
      maximum_number_of_rows=-1,
      maximum_number_of_columns=-1)
  else
    println("  (no entries found)")
  end

  # ─── Summary ─────────────────────────────────────────────────────────
  println("\n" * "=" ^ 100)
  println(" Summary & Benchmarking")
  println("=" ^ 100)

  n_total = length(results)
  n_ok = count(r -> r.hodge_match && r.inv_match, results)
  n_mismatch = count(
    r -> r.error_msg === nothing && !(r.hodge_match && r.inv_match), results
  )
  n_error = count(r -> r.error_msg !== nothing, results)
  n_undetermined = count(r -> r.underdetermined, results)

  println("  Entries processed:          $n_total")
  println("  Matches (correct):          $n_ok")
  println("  Mismatches:                 $n_mismatch")
  println("  Computation errors:         $n_error")
  println("  Underdetermined (symbolic): $n_undetermined")
  println()
  @printf("  Wall-clock time: %.1f seconds\n", elapsed_total)
  per_entry_times = [r.elapsed for r in results if r.error_msg === nothing]
  if !isempty(per_entry_times)
    @printf("  Sum of per-entry times:     %.1f seconds\n", sum(per_entry_times))
    @printf(
      "  Average time per entry:     %.3f seconds\n",
      sum(per_entry_times) / length(per_entry_times)
    )
    @printf(
      "  Time range (per entry):     %.3f – %.2f seconds\n",
      minimum(per_entry_times),
      maximum(per_entry_times)
    )
  end

  # ─── Underdetermined entries ─────────────────────────────────────────
  if n_undetermined > 0
    println("\n─── Underdetermined entries (symbolic Hodge numbers) ───")
    for r in results
      r.underdetermined || continue
      H = r.hodge_computed !== nothing ? r.hodge_computed : r.hodge_reference
      vars = Set{Int}()
      for e in H
        union!(vars, keys(e.coeffs))
      end
      var_names = join(["z_$v" for v in sort!(collect(vars))], ", ")
      @printf("  %-7s  %-25s  free variables: %s\n", r.id, r.ambient_desc, var_names)
    end
  end

  # ─── Mismatches ──────────────────────────────────────────────────────
  if n_mismatch > 0
    println("\n─── Mismatches ───")
    for r in results
      r.error_msg === nothing && !(r.hodge_match && r.inv_match) || continue
      println("  $(r.id):")
      if !r.hodge_match
        println("    Hodge: MISMATCH")
        println("    Computed:  $(r.hodge_computed)")
        println("    Reference: $(r.hodge_reference)")
      end
      if !r.inv_match
        @printf("    h⁰(-K): computed=%s  reference=%d\n",
          r.h0_antiK === nothing ? "?" : string(r.h0_antiK), r.ref_h0_antiK)
        @printf("    (-K)⁴:  computed=%s  reference=%d\n",
          r.antiK4 === nothing ? "?" : string(r.antiK4), r.ref_antiK4)
        @printf("    χ(T):   computed=%s  reference=%d\n",
          r.chi_T === nothing ? "?" : string(r.chi_T), r.ref_chi_T)
      end
    end
  end

  # ─── Errors ──────────────────────────────────────────────────────────
  if n_error > 0
    println("\n─── Errors ───")
    for r in results
      r.error_msg === nothing && continue
      @printf("  %-7s: %s\n", r.id, r.error_msg)
    end
  end

  if n_mismatch == 0 && n_error == 0
    println("\n  ✅ All $(n_total) entries matched the reference data successfully.")
    println("     Reference: https://zenodo.org/doi/10.5281/zenodo.17742466")
  end

  results
end

# ═══════════════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════════════

if abspath(PROGRAM_FILE) == @__FILE__
  function _cli_main()
    ids = nothing

    i = 1
    while i <= length(ARGS)
      if ARGS[i] == "--ids" && i < length(ARGS)
        ids = String.(split(ARGS[i + 1], ","))
        i += 2
      else
        i += 1
      end
    end

    main(; ids=ids)
  end

  _cli_main()
end
