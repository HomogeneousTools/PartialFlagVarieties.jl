# ═══════════════════════════════════════════════════════════════════════════════
#  FanoFourfoldsK3.jl — FK3 fourfolds from arXiv:2111.13030
#
#  Computes Hodge numbers for the 64 Fano fourfolds of K3 type from
#    Bernardara–Fatighenti–Manivel–Tanturri,
#    "Fano fourfolds of K3 type" (arXiv:2111.13030).
#
#  The descriptions are extracted at runtime from the family headings in the
#  paper. Ambient/bundle data are located inside examples/FanoFourfolds.json by
#  matching the Appendix C invariants together with the ambient factorization.
#
#  The script then recomputes the Hodge diamond with PartialFlagVarieties and
#  reports:
#    • matches with the paper table,
#    • numerical discrepancies,
#    • indeterminacies left by `hodge_numbers_symbolic`.
#
#  Usage:
#    julia --project=. examples/FanoFourfoldsK3.jl
#    julia --project=. examples/FanoFourfoldsK3.jl K3-29 K3-30 R-64
#
#  Last verified comparison status:
#    • 58 matches
#    • 5 symbolic indeterminacies: C-2, C-11, R-61, R-63, R-64
#    • 1 determined discrepancy: K3-37
#
#  K3-37 deserves special attention: unlike the indeterminate families, its
#  Hodge numbers are fully determined by the current code. The script computes
#    (h11, h21, h22) = (3, 2, 22)
#  for
#    𝒵(ℙ1×ℙ1×ℙ5, 𝒪(0,0,2) ⊕ 𝒪(0,0,2) ⊕ 𝒪(1,1,1)),
#  while Appendix C lists
#    (h11, h21, h22) = (3, 2, 33).
#  This therefore looks like a genuine discrepancy rather than a symbolic gap.
# ═══════════════════════════════════════════════════════════════════════════════

using Downloads
using JSON
using Lie: rank
using PartialFlagVarieties
using PrettyTables
using Printf

const PAPER_HTML_URL = "https://arxiv.org/html/2111.13030v2"
const LOCAL_DB_PATH = joinpath(@__DIR__, "FanoFourfolds.json")

normalize_label(s::AbstractString) = replace(strip(s), '–' => '-', '—' => '-', '−' => '-')

Base.@kwdef struct FamilySpec
  label::String
  h11::Int
  h21::Int
  h22::Int
  h0antiK::Int
  volume::Int
  minus_chiT::Int
end

Base.@kwdef struct FamilyResult
  label::String
  description::String
  ambient::String
  computed_h11
  computed_h21
  computed_h22
  expected_h11::Int
  expected_h21::Int
  expected_h22::Int
  status::Symbol
  note::String
end

const FAMILY_SPECS = let
  raw = """
  C-1|2|0|28|36|144|28
  C-2|2|1|23|28|99|29
  C-3|2|0|22|45|192|19
  C-4|2|0|23|40|163|24
  C-5|2|0|23|39|161|21
  C-6|2|0|21|49|211|19
  C-7|2|0|26|24|86|24
  C-8|2|0|27|30|115|26
  C-9|2|0|23|27|101|21
  C-10|3|0|38|20|63|28
  C-11|3|1|30|24|81|31
  C-12|3|0|30|27|99|27
  C-13|3|0|28|34|134|25
  C-14|3|0|30|30|113|28
  C-15|3|0|23|35|141|18
  C-16|3|0|23|42|176|18
  C-17|3|0|24|37|149|21
  GM-18|2|0|30|22|74|30
  GM-19|2|0|27|23|80|26
  GM-20|2|0|28|26|94|28
  GM-21|2|0|24|30|114|24
  GM-22|2|0|23|33|130|22
  GM-23|2|0|24|30|115|23
  K3-24|2|0|32|19|60|31
  K3-25|2|0|28|21|70|27
  K3-26|2|0|24|23|80|23
  K3-27|2|0|28|32|124|28
  K3-28|2|0|22|25|90|21
  K3-29|2|0|22|39|160|21
  K3-30|2|0|22|39|160|21
  K3-31|2|0|28|23|80|27
  K3-32|2|0|22|29|110|21
  K3-33|2|0|24|27|100|23
  K3-34|2|0|23|26|95|22
  K3-35|2|0|24|25|90|23
  K3-36|3|5|22|19|60|23
  K3-37|3|2|33|23|80|20
  K3-38|3|0|22|27|100|18
  K3-39|3|0|24|29|110|20
  K3-40|3|0|32|22|74|29
  K3-41|3|0|22|39|160|18
  K3-42|3|0|27|28|104|24
  K3-43|3|0|24|32|125|20
  K3-44|3|0|30|25|89|27
  K3-45|3|0|24|31|119|21
  K3-46|3|0|24|33|130|20
  K3-47|3|0|26|27|100|22
  K3-48|3|0|22|31|120|18
  K3-49|3|0|24|27|100|20
  K3-50|3|0|22|43|180|18
  K3-51|3|0|30|23|80|26
  K3-52|3|0|24|29|110|20
  K3-53|3|0|23|36|145|19
  K3-54|3|0|25|34|133|23
  K3-55|3|0|22|37|150|18
  K3-56|4|0|24|35|140|17
  K3-57|4|0|24|34|134|18
  K3-58|4|0|28|28|104|22
  K3-59|4|0|24|31|120|17
  K3-60|5|0|26|31|120|16
  R-61|2|0|25|27|99|25
  R-62|2|0|23|33|129|23
  R-63|3|2|22|31|120|20
  R-64|3|5|22|25|90|23
  """

  map(filter(line -> !isempty(strip(line)), split(strip(raw), '\n'))) do line
    fields = split(strip(line), '|')
    FamilySpec(
      label = fields[1],
      h11 = parse(Int, fields[2]),
      h21 = parse(Int, fields[3]),
      h22 = parse(Int, fields[4]),
      h0antiK = parse(Int, fields[5]),
      volume = parse(Int, fields[6]),
      minus_chiT = parse(Int, fields[7]),
    )
  end
end

const SELECTED_LABELS = Set{String}(normalize_label(arg) for arg in ARGS)

function remove_invisibles(s::AbstractString)
  bad = Set(['\u2062', '\u2061', '\u2060', '\ufeff'])
  String([c for c in s if !(c in bad)])
end

function normalize_factor_text(s::AbstractString)
  sub_map = Dict(
    '₀' => '0', '₁' => '1', '₂' => '2', '₃' => '3', '₄' => '4',
    '₅' => '5', '₆' => '6', '₇' => '7', '₈' => '8', '₉' => '9',
    '⁰' => '0', '¹' => '1', '²' => '2', '³' => '3', '⁴' => '4',
    '⁵' => '5', '⁶' => '6', '⁷' => '7', '⁸' => '8', '⁹' => '9'
  )
  t = remove_invisibles(strip(s))
  t = join(get(sub_map, c, c) for c in t)
  replace(t, r"\s+" => "")
end

function top_level_split(s::AbstractString, delim::Char)
  parts = String[]
  depth = 0
  start = firstindex(s)
  i = firstindex(s)
  while i <= lastindex(s)
    c = s[i]
    if c == '('
      depth += 1
    elseif c == ')'
      depth -= 1
    elseif c == delim && depth == 0
      push!(parts, strip(s[start:prevind(s, i)]))
      start = nextind(s, i)
    end
    i = nextind(s, i)
  end
  push!(parts, strip(s[start:end]))
  parts
end

function extract_formula_line(block::AbstractString)
  start = findfirst('𝒵', block)
  start === nothing && return ""
  stop = findnext('.', block, start)
  formula = stop === nothing ? strip(block[start:end]) : strip(block[start:prevind(block, stop)])
  second_z = findnext('𝒵', formula, nextind(formula, firstindex(formula)))
  second_z === nothing || (formula = strip(formula[firstindex(formula):prevind(formula, second_z)]))
  for marker in (" \\math", " script_", " \\mathscr", " \\displaystyle")
    pos = findfirst(marker, formula)
    pos === nothing || (formula = strip(formula[firstindex(formula):prevind(formula, first(pos))]))
  end
  formula = replace(formula, r"\s+" => " ")
  formula = replace(formula, r"\s*\(\s*" => "(")
  formula = replace(formula, r"\s*\)\s*" => ")")
  formula = replace(formula, r"\s*,\s*" => ",")
  formula = replace(formula, r"\s*⊕\s*" => " ⊕ ")
  formula = replace(formula, r"\s*×\s*" => "×")
  formula
end

function fetch_paper_descriptions()
  html_path = Downloads.download(PAPER_HTML_URL)
  html = read(html_path, String)
  text = replace(html, r"<script[\s\S]*?</script>" => " ")
  text = replace(text, r"<style[\s\S]*?</style>" => " ")
  text = replace(text, r"<[^>]+>" => " ")
  text = replace(text, "&nbsp;" => " ", "&amp;" => "&")
  text = remove_invisibles(text)
  text = replace(text, r"\s+" => " ")

  descs = Dict{String,String}()
  for spec in FAMILY_SPECS
    label_dash = replace(spec.label, "-" => "–")
    pat = Regex("Fano\\s+" * escape_string(label_dash) * "\\s*\\.")
    m = match(pat, text)
    m === nothing && continue
    stop = m.offset
    for _ in 1:1200
      stop == lastindex(text) && break
      stop = nextind(text, stop)
    end
    tail = text[m.offset:stop]
    formula = extract_formula_line(tail)
    !isempty(formula) && (descs[spec.label] = formula)
  end
  descs
end

function extract_ambient_part(description::AbstractString)
  start = findfirst('(', description)
  start === nothing && error("Cannot find '(' in description: $description")
  depth = 0
  i = nextind(description, start)
  comma_at = nothing
  while i <= lastindex(description)
    c = description[i]
    if c == '('
      depth += 1
    elseif c == ')'
      depth -= 1
    elseif c == ',' && depth == 0
      comma_at = i
      break
    end
    i = nextind(description, i)
  end
  comma_at === nothing && error("Cannot isolate ambient part in: $description")
  strip(description[nextind(description, start):prevind(description, comma_at)])
end

function parse_projective_power(factor::AbstractString)
  m = match(r"^\(ℙ1\)(\d+)$", factor)
  m === nothing && return nothing
  n = parse(Int, m.captures[1])
  [[1, 2] for _ in 1:n]
end

function parse_ambient(description::AbstractString)
  ambient_part = extract_ambient_part(description)
  factors = top_level_split(ambient_part, '×')
  parsed = Vector{Vector{Int}}()

  for factor0 in factors
    factor = normalize_factor_text(factor0)
    repeated = parse_projective_power(factor)
    if repeated !== nothing
      append!(parsed, repeated)
      continue
    end

    factor = replace(factor, r"_[0-9]+$" => "")
    factor = replace(factor, r"ℙ([0-9]+)[0-9]+$" => s"ℙ\1")
    factor = replace(factor, r"Gr\(([^\)]*)\)[0-9]+$" => s"Gr(\1)")

    if startswith(factor, 'ℙ')
      m = match(r"^ℙ(\d+)$", factor)
      m === nothing && error("Unsupported projective factor: $factor0")
      digits = m.captures[1]
      if occursin(r"[₀₁₂₃₄₅₆₇₈₉]", factor0) && length(digits) > 1
        digits = string(first(digits))
      end
      push!(parsed, [1, parse(Int, digits) + 1])
    elseif startswith(factor, "Gr(")
      m = match(r"^Gr\((\d+),(\d+)\)$", factor)
      m === nothing && error("Unsupported Grassmannian factor: $factor0")
      push!(parsed, [parse(Int, m.captures[1]), parse(Int, m.captures[2])])
    elseif startswith(factor, "Fl(")
      m = match(r"^Fl\((\d+(?:,\d+)*)\)$", factor)
      m === nothing && error("Unsupported flag factor: $factor0")
      nums = parse.(Int, split(m.captures[1], ','))
      push!(parsed, nums)
    else
      error("Unsupported ambient factor: $factor0")
    end
  end

  parsed
end

function gl_weight_to_omega_flag(ks::Vector{Int}, n::Int, w::Vector{Int})
  eps = Int[]
  append!(eps, w[n - ks[1] + 1:n])
  for i in 2:length(ks)
    append!(eps, w[n - ks[i] + 1:n - ks[i - 1]])
  end
  append!(eps, w[1:n - ks[end]])
  [eps[i] - eps[i + 1] for i in 1:n - 1]
end

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

  partial_flag_variety(DT, Tuple(all_marks))
end

function build_bundle(X, factors::Vector, bundle_data::Vector)
  mdt = marked_dynkin_type(X)

  factor_ks = Vector{Int}[]
  factor_ns = Int[]
  for f in factors
    fv = Int.(f)
    push!(factor_ks, fv[1:end - 1])
    push!(factor_ns, fv[end])
  end

  summands = IrrepLevi[]
  isempty(bundle_data) && return zero_bundle(X)
  ds = bundle_data[1]
  (ds isa Vector && isempty(ds)) && return zero_bundle(X)

  for summand in ds
    omega_coords = Int[]
    for (j, factor_weight) in enumerate(summand)
      w = Int.(factor_weight[1])
      omega = gl_weight_to_omega_flag(factor_ks[j], factor_ns[j], w)
      append!(omega_coords, omega)
    end
    lam = WeightLatticeElem(dynkin_type(X), omega_coords)
    push!(summands, IrrepLevi(mdt, lam))
  end

  CompletelyReducibleBundle(X, summands)
end

function build_bundle_from_factor_weights(X, factors::Vector, summands_weights)
  mdt = marked_dynkin_type(X)

  factor_ks = Vector{Int}[]
  factor_ns = Int[]
  for f in factors
    fv = Int.(f)
    push!(factor_ks, fv[1:end - 1])
    push!(factor_ns, fv[end])
  end

  summands = IrrepLevi[]
  for summand in summands_weights
    omega_coords = Int[]
    for (j, w_any) in enumerate(summand)
      w = Int.(w_any)
      omega = gl_weight_to_omega_flag(factor_ks[j], factor_ns[j], w)
      append!(omega_coords, omega)
    end
    lam = WeightLatticeElem(dynkin_type(X), omega_coords)
    push!(summands, IrrepLevi(mdt, lam))
  end

  CompletelyReducibleBundle(X, summands)
end

projective_O_weight(dim::Int, degree::Int) = vcat(zeros(Int, dim), [degree])
projective_Q_weight(dim::Int) = vcat(ones(Int, dim - 1), [0, 1])
grassmannian_Udual_weight(k::Int, n::Int) = begin
  w = zeros(Int, n)
  w[n - k + 1] = 1
  w
end
grassmannian_O1_weight(k::Int, n::Int) = vcat(zeros(Int, n - k), ones(Int, k))

const MANUAL_OVERRIDES = Dict{String,NamedTuple}(
  "C-5" => (
    ambient = Any[[1, 6], [2, 8]],
    summands = Any[
      Any[projective_O_weight(5, 1), grassmannian_Udual_weight(2, 8)],
      Any[projective_O_weight(5, 1), grassmannian_O1_weight(2, 8)],
      Any[projective_Q_weight(5), grassmannian_Udual_weight(2, 8)],
    ],
  ),
  "C-10" => (
    ambient = Any[[1, 2], [1, 2], [1, 6]],
    summands = Any[
      Any[projective_O_weight(1, 0), projective_O_weight(1, 0), projective_O_weight(5, 3)],
      Any[projective_O_weight(1, 0), projective_O_weight(1, 1), projective_O_weight(5, 1)],
      Any[projective_O_weight(1, 1), projective_O_weight(1, 0), projective_O_weight(5, 1)],
    ],
  ),
  "C-15" => (
    ambient = Any[[1, 3], [1, 3], [1, 6]],
    summands = Any[
      Any[projective_Q_weight(2), projective_O_weight(2, 0), projective_O_weight(5, 1)],
      Any[projective_O_weight(2, 0), projective_Q_weight(2), projective_O_weight(5, 1)],
      Any[projective_O_weight(2, 1), projective_O_weight(2, 1), projective_O_weight(5, 1)],
    ],
  ),
  "K3-37" => (
    ambient = Any[[1, 2], [1, 2], [1, 6]],
    summands = Any[
      Any[projective_O_weight(1, 0), projective_O_weight(1, 0), projective_O_weight(5, 2)],
      Any[projective_O_weight(1, 0), projective_O_weight(1, 0), projective_O_weight(5, 2)],
      Any[projective_O_weight(1, 1), projective_O_weight(1, 1), projective_O_weight(5, 1)],
    ],
  ),
  "K3-40" => (
    ambient = Any[[1, 2], [1, 2], [1, 6]],
    summands = Any[
      Any[projective_O_weight(1, 0), projective_O_weight(1, 0), projective_O_weight(5, 2)],
      Any[projective_O_weight(1, 0), projective_O_weight(1, 1), projective_O_weight(5, 1)],
      Any[projective_O_weight(1, 1), projective_O_weight(1, 0), projective_O_weight(5, 2)],
    ],
  ),
  "K3-43" => (
    ambient = Any[[1, 3], [1, 3], [1, 5]],
    summands = Any[
      Any[projective_O_weight(2, 1), projective_O_weight(2, 0), projective_O_weight(4, 1)],
      Any[projective_O_weight(2, 1), projective_O_weight(2, 1), projective_O_weight(4, 1)],
      Any[projective_O_weight(2, 0), projective_Q_weight(2), projective_O_weight(4, 1)],
    ],
  ),
  "K3-45" => (
    ambient = Any[[1, 3], [1, 3], [1, 5]],
    summands = Any[
      Any[projective_O_weight(2, 0), projective_O_weight(2, 1), projective_O_weight(4, 2)],
      Any[projective_O_weight(2, 1), projective_O_weight(2, 1), projective_O_weight(4, 0)],
      Any[projective_Q_weight(2), projective_O_weight(2, 0), projective_O_weight(4, 1)],
    ],
  ),
  "K3-50" => (
    ambient = Any[[1, 4], [1, 5], [1, 5]],
    summands = Any[
      Any[projective_O_weight(3, 1), projective_O_weight(4, 1), projective_O_weight(4, 1)],
      Any[projective_Q_weight(3), projective_O_weight(4, 1), projective_O_weight(4, 0)],
      Any[projective_Q_weight(3), projective_O_weight(4, 0), projective_O_weight(4, 1)],
    ],
  ),
  "K3-58" => (
    ambient = Any[[1, 2], [1, 2], [1, 2], [1, 4]],
    summands = Any[
      Any[projective_O_weight(1, 0), projective_O_weight(1, 0), projective_O_weight(1, 1), projective_O_weight(3, 1)],
      Any[projective_O_weight(1, 1), projective_O_weight(1, 1), projective_O_weight(1, 0), projective_O_weight(3, 2)],
    ],
  ),
  "K3-59" => (
    ambient = Any[[1, 2], [1, 2], [1, 3], [1, 3]],
    summands = Any[
      Any[projective_O_weight(1, 0), projective_O_weight(1, 0), projective_O_weight(2, 1), projective_O_weight(2, 1)],
      Any[projective_O_weight(1, 1), projective_O_weight(1, 1), projective_O_weight(2, 1), projective_O_weight(2, 1)],
    ],
  ),
)

const PROBLEMATIC_FAMILY_SUMMARIES = Dict(
  "C-2" =>
    "Symbolic case: the ambient/model match is clear, but `hodge_numbers_symbolic` leaves h21 = x_12 and h22 = 21 + 2*x_12 unresolved.",
  "C-11" =>
    "Symbolic case: the long-exact-sequence stage does not eliminate all affine variables, so both h21 and h22 remain undetermined.",
  "K3-37" =>
    "Determined discrepancy: this family is modeled manually as 𝒵(ℙ1×ℙ1×ℙ5, 𝒪(0,0,2) ⊕ 𝒪(0,0,2) ⊕ 𝒪(1,1,1)), and the code consistently returns (3,2,22), not the Appendix C value (3,2,33).",
  "R-61" =>
    "Symbolic case: even h11 is not forced numerically; the output depends on the free combination x_13 - x_9.",
  "R-63" =>
    "Symbolic case: the unresolved affine parameters affect all three reported entries, so the expected Picard-rank jump cannot be confirmed inside the current code.",
  "R-64" =>
    "Symbolic case: same pattern as R-63, with free variables surviving in h11, h21, and h22.",
)

function maybe_int(x)
  x isa Integer && return Int(x)
  x isa AbstractString || return nothing
  s = strip(x)
  occursin(r"^-?\d+$", s) || return nothing
  parse(Int, s)
end

entry_h11(entry) = maybe_int(entry["hodge"][2][2])
entry_h21(entry) = maybe_int(entry["hodge"][2][3])
entry_h22(entry) = maybe_int(entry["hodge"][3][3])
entry_minus_chiT(entry) = -Int(entry["tangent"])

function ambient_to_string(ambient)
  join(map(ambient) do factor
    fv = Int.(factor)
    if length(fv) == 2 && fv[1] == 1
      "P$(fv[2] - 1)"
    elseif length(fv) == 2
      "Gr($(fv[1]),$(fv[2]))"
    else
      "Fl(" * join(string.(fv), ",") * ")"
    end
  end, " × ")
end

function affine_expr_string(e::AffineExpr)
  string(e)
end

function matrix_entry_exprs(H::AbstractMatrix{AffineExpr})
  (
    h11 = H[2, 2],
    h21 = H[3, 2],
    h22 = H[3, 3],
  )
end

function all_relevant_determined(exprs)
  is_determined(exprs.h11) && is_determined(exprs.h21) && is_determined(exprs.h22)
end

function format_expr(e::AffineExpr)
  is_determined(e) ? string(e.constant) : affine_expr_string(e)
end

function select_specs()
  isempty(SELECTED_LABELS) && return FAMILY_SPECS
  selected = [spec for spec in FAMILY_SPECS if spec.label in SELECTED_LABELS]
  isempty(selected) && error("No labels matched ARGS = $(join(ARGS, ", "))")
  selected
end

function load_database()
  JSON.parsefile(LOCAL_DB_PATH)
end

function find_candidates(spec::FamilySpec, description::String, db)
  ambient = parse_ambient(description)
  strict = [entry for entry in db if
    entry["ambient"] == ambient &&
    Int(entry["anticanonical"]) == spec.h0antiK &&
    Int(entry["volume"]) == spec.volume &&
    entry_minus_chiT(entry) == spec.minus_chiT &&
    (entry_h11(entry) === nothing || entry_h11(entry) == spec.h11) &&
    (entry_h21(entry) === nothing || entry_h21(entry) == spec.h21) &&
    (entry_h22(entry) === nothing || entry_h22(entry) == spec.h22)
  ]
  !isempty(strict) && return strict

  [entry for entry in db if
    entry["ambient"] == ambient &&
    Int(entry["anticanonical"]) == spec.h0antiK &&
    Int(entry["volume"]) == spec.volume &&
    entry_minus_chiT(entry) == spec.minus_chiT
  ]
end

function compute_family(spec::FamilySpec, description::String, db)
  candidates = find_candidates(spec, description, db)
  ambient_str = ambient_to_string(parse_ambient(description))

  if isempty(candidates) && haskey(MANUAL_OVERRIDES, spec.label)
    override = MANUAL_OVERRIDES[spec.label]
    X = build_ambient(override.ambient)
    E = build_bundle_from_factor_weights(X, override.ambient, override.summands)
    ambient_str = ambient_to_string(override.ambient)
    note_prefix = "manual override"
  elseif isempty(candidates)
    return FamilyResult(
      label = spec.label,
      description = description,
      ambient = ambient_str,
      computed_h11 = "—",
      computed_h21 = "—",
      computed_h22 = "—",
      expected_h11 = spec.h11,
      expected_h21 = spec.h21,
      expected_h22 = spec.h22,
      status = :unmatched,
      note = "No matching entry found in examples/FanoFourfolds.json",
    )
  elseif length(candidates) > 1
    return FamilyResult(
      label = spec.label,
      description = description,
      ambient = ambient_str,
      computed_h11 = "—",
      computed_h21 = "—",
      computed_h22 = "—",
      expected_h11 = spec.h11,
      expected_h21 = spec.h21,
      expected_h22 = spec.h22,
      status = :ambiguous,
      note = "Matched $(length(candidates)) database entries with the same invariants",
    )
  else
    entry = only(candidates)
    X = build_ambient(entry["ambient"])
    E = build_bundle(X, entry["ambient"], entry["bundle"])
    note_prefix = ""
  end

  Z = zero_locus(E)
  Hsym = hodge_numbers_symbolic(Z)
  exprs = matrix_entry_exprs(Hsym)

  if all_relevant_determined(exprs)
    got_h11 = Int(exprs.h11.constant)
    got_h21 = Int(exprs.h21.constant)
    got_h22 = Int(exprs.h22.constant)
    ok = got_h11 == spec.h11 && got_h21 == spec.h21 && got_h22 == spec.h22
    return FamilyResult(
      label = spec.label,
      description = description,
      ambient = ambient_str,
      computed_h11 = got_h11,
      computed_h21 = got_h21,
      computed_h22 = got_h22,
      expected_h11 = spec.h11,
      expected_h21 = spec.h21,
      expected_h22 = spec.h22,
      status = ok ? :match : :discrepancy,
      note = ok ? note_prefix : "Computed values disagree with Appendix C",
    )
  end

  FamilyResult(
    label = spec.label,
    description = description,
    ambient = ambient_str,
    computed_h11 = format_expr(exprs.h11),
    computed_h21 = format_expr(exprs.h21),
    computed_h22 = format_expr(exprs.h22),
    expected_h11 = spec.h11,
    expected_h21 = spec.h21,
    expected_h22 = spec.h22,
    status = :indeterminate,
    note = isempty(note_prefix) ? "Symbolic ambiguity in hodge_numbers_symbolic" : note_prefix * "; symbolic ambiguity in hodge_numbers_symbolic",
  )
end

function show_table(results)
  rows = map(results) do r
    [
      r.label,
      r.ambient,
      r.computed_h11,
      r.expected_h11,
      r.computed_h21,
      r.expected_h21,
      r.computed_h22,
      r.expected_h22,
      String(r.status),
      isempty(r.note) ? "—" : r.note,
    ]
  end
  data = permutedims(hcat(rows...), (2, 1))
  pretty_table(
    data;
    column_labels = ["label", "ambient", "h11", "exp h11", "h21", "exp h21", "h22", "exp h22", "status", "note"],
    alignment = [:l, :l, :r, :r, :r, :r, :r, :r, :l, :l],
    fit_table_in_display_horizontally = false,
    fit_table_in_display_vertically = false,
  )
end

function show_issues(results)
  bad = [r for r in results if r.status != :match]
  println()
  if isempty(bad)
    println("No discrepancies or indeterminacies found.")
    return
  end

  println("Discrepancies / indeterminacies:")
  for r in bad
    println("- $(r.label): $(isempty(r.note) ? String(r.status) : r.note)")
    println("  $(r.description)")
    println("  computed = (h11=$(r.computed_h11), h21=$(r.computed_h21), h22=$(r.computed_h22))")
    println("  expected = (h11=$(r.expected_h11), h21=$(r.expected_h21), h22=$(r.expected_h22))")
    if haskey(PROBLEMATIC_FAMILY_SUMMARIES, r.label)
      println("  summary  = $(PROBLEMATIC_FAMILY_SUMMARIES[r.label])")
    end
  end
end

function main()
  specs = select_specs()
  db = load_database()
  descriptions = fetch_paper_descriptions()

  missing = [spec.label for spec in specs if !haskey(descriptions, spec.label)]
  !isempty(missing) && error("Failed to extract descriptions for: $(join(missing, ", "))")

  println("Computing $(length(specs)) FK3 family/families from arXiv:2111.13030 ...")
  println()

  results = FamilyResult[]
  for (i, spec) in enumerate(specs)
    desc = descriptions[spec.label]
    print(@sprintf("[%2d/%2d] %-6s ", i, length(specs), spec.label))
    flush(stdout)
    t0 = time()
    res = compute_family(spec, desc, db)
    push!(results, res)
    println(@sprintf("%s (%.2fs)", String(res.status), time() - t0))
  end

  println()
  show_table(results)
  show_issues(results)
end

main()