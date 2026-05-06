# ═══════════════════════════════════════════════════════════════════════════════
#  Labels.jl — ZeroLocus62 label encoding/decoding for G/P and zero loci
#
#  Provides bidirectional conversion between PartialFlagVariety / ZeroLocus
#  objects and ZeroLocus62 canonical label strings.
# ═══════════════════════════════════════════════════════════════════════════════

export zerolocus62_label

using ZeroLocus62: Factor, encode_label, decode_label

# ═══════════════════════════════════════════════════════════════════════════════
#  Internal helpers: DynkinType ↔ Factor conversion
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _simple_dynkin_family(DT::Type{<:SimpleDynkinType}) -> Char

Return the Dynkin family letter for a simple Dynkin type.
"""
function _simple_dynkin_family(::Type{DT}) where {DT<:SimpleDynkinType}
  DT <: TypeA && return 'A'
  DT <: TypeB && return 'B'
  DT <: TypeC && return 'C'
  DT <: TypeD && return 'D'
  DT <: TypeE && return 'E'
  DT == TypeF4 && return 'F'
  DT == TypeG2 && return 'G'
  error("Unknown simple Dynkin type: $DT")
end

"""
    _flatten_simple_factors(DT::Type{<:DynkinType}) -> Vector{DataType}

Recursively flatten a (possibly nested) `ProductDynkinType` into an ordered
list of simple factors.
"""
function _flatten_simple_factors(::Type{DT}) where {DT<:SimpleDynkinType}
  DataType[DT]
end

function _flatten_simple_factors(::Type{DT}) where {DT<:ProductDynkinType}
  result = DataType[]
  for T in DT.parameters[1].parameters
    append!(result, _flatten_simple_factors(T))
  end
  result
end

"""
    _mdt_to_factors(mdt::MarkedDynkinType) -> Vector{Factor}

Convert a `MarkedDynkinType` to a vector of `ZeroLocus62.Factor` objects, one
per irreducible Dynkin factor. For a simple type this is a single factor; for a
product type the marked nodes are partitioned by cumulative rank offsets.
"""
function _mdt_to_factors(mdt::MarkedDynkinType)
  DT = dynkin_type(mdt)
  simple_factors = _flatten_simple_factors(DT)
  marked = marked_nodes(mdt)

  factors = Factor[]
  offset = 0
  for sf in simple_factors
    r = rank(sf)
    local_marked = [m - offset for m in marked if offset < m <= offset + r]
    mask = foldl(|, (1 << (m - 1) for m in local_marked); init=0)
    mask == 0 && throw(
      ArgumentError(
        "Factor $(_simple_dynkin_family(sf))$r has no marked nodes; " *
        "ZeroLocus62 requires at least one marked node per factor.",
      ),
    )
    push!(factors, Factor(_simple_dynkin_family(sf), r, mask))
    offset += r
  end

  factors
end

"""
    _char_to_symbol(c::Char) -> Symbol

Convert a Dynkin family character to the corresponding Symbol for
`_symbol_to_simple_type`.
"""
_char_to_symbol(c::Char) = Symbol(c)

"""
    _factors_to_mdt(factors::Vector{Factor}) -> MarkedDynkinType

Convert a vector of `ZeroLocus62.Factor` objects back to a `MarkedDynkinType`.
"""
function _factors_to_mdt(factors::Vector{Factor})
  length(factors) >= 1 || throw(ArgumentError("Need at least one factor"))

  simple_types = DataType[
    _symbol_to_simple_type(_char_to_symbol(f.group), f.rank) for f in factors
  ]

  if length(simple_types) == 1
    DT = simple_types[1]
  else
    DT = simple_types[1]
    for i in 2:length(simple_types)
      DT = ProductDynkinType{Tuple{DT,simple_types[i]}}
    end
  end

  global_marked = Int[]
  offset = 0
  for f in factors
    for node in 1:f.rank
      if ((f.mask >> (node - 1)) & 1) == 1
        push!(global_marked, node + offset)
      end
    end
    offset += f.rank
  end

  MarkedDynkinType(DT, Tuple(sort(global_marked)))
end

"""
    _weight_to_summand_row(λ::WeightLatticeElem, factor_ranks::Vector{Int}) -> Vector{Vector{Int}}

Split a weight's coefficient vector into per-factor weight vectors.
"""
function _weight_to_summand_row(λ::WeightLatticeElem, factor_ranks::Vector{Int})
  coeffs = Int[c for c in Lie.coefficients(λ)]
  row = Vector{Vector{Int}}()
  offset = 1
  for r in factor_ranks
    push!(row, coeffs[offset:(offset + r - 1)])
    offset += r
  end
  row
end

"""
    _summand_row_to_weight(row::Vector{Vector{Int}}, DT::DataType) -> WeightLatticeElem

Concatenate per-factor weight vectors into a single `WeightLatticeElem`.
"""
function _summand_row_to_weight(row::Vector{Vector{Int}}, DT::DataType)
  coeffs = reduce(vcat, row)
  WeightLatticeElem(DT, coeffs)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Export: PartialFlagVariety / ZeroLocus → label
# ═══════════════════════════════════════════════════════════════════════════════

"""
    zerolocus62_label(X::PartialFlagVariety) -> String

Encode a partial flag variety as a ZeroLocus62 label (ambient-only, no bundle).

This is the label-level counterpart of the ordinary constructor
[`PartialFlagVariety(label::AbstractString)`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> zerolocus62_label(projective_space(1))
"1"

julia> zerolocus62_label(Gr(2, 4))
"31"

julia> zerolocus62_label(projective_space(3))
"30"
```
"""
function zerolocus62_label(X::PartialFlagVariety)
  factors = _mdt_to_factors(marked_dynkin_type(X))
  encode_label(factors, Vector{Vector{Vector{Int}}}())
end

"""
    zerolocus62_label(Z::ZeroLocus) -> String

Encode a zero locus (ambient variety + defining bundle) as a ZeroLocus62 label.

Each irreducible summand of the defining bundle contributes one summand row
whose entries are the fundamental-weight coefficients of its ``P``-dominant weight,
split by ambient factor.

The label records the ambient variety and defining bundle only. It does not
encode any proof that the corresponding bundle admits a regular section.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(1);

julia> zerolocus62_label(zero_locus(line_bundle(X, 1)))
"1.0"
```
"""
function zerolocus62_label(Z::ZeroLocus)
  X = ambient_variety(Z)
  factors = _mdt_to_factors(marked_dynkin_type(X))
  factor_ranks = [f.rank for f in factors]

  E = defining_bundle(Z)
  summands = Vector{Vector{Vector{Int}}}()
  for comp in components(E)
    λ = p_dominant_weight(comp)
    push!(summands, _weight_to_summand_row(λ, factor_ranks))
  end

  encode_label(factors, summands)
end

"""
    zerolocus62_label(E::CompletelyReducibleBundle) -> String

Encode an equivariant bundle on a partial flag variety as a ZeroLocus62 label
(ambient + bundle summands). This does not require the bundle to define a
valid zero locus.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(1);

julia> zerolocus62_label(direct_sum(structure_sheaf(X), line_bundle(X, 1)))
"1.x10"
```
"""
function zerolocus62_label(E::CompletelyReducibleBundle)
  X = variety(E)
  factors = _mdt_to_factors(marked_dynkin_type(X))
  factor_ranks = [f.rank for f in factors]

  summands = Vector{Vector{Vector{Int}}}()
  for comp in components(E)
    λ = p_dominant_weight(comp)
    push!(summands, _weight_to_summand_row(λ, factor_ranks))
  end

  encode_label(factors, summands)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Import: label → PartialFlagVariety / ZeroLocus
# ═══════════════════════════════════════════════════════════════════════════════

"""
    PartialFlagVariety(label::AbstractString) -> PartialFlagVariety

Construct a partial flag variety from a ZeroLocus62 label. Any bundle part
in the label is ignored.

This constructor is intentionally different from
[`partial_flag_variety(s, marked)`](@ref): there, the string `s` is a Dynkin
type such as `"A4"`, whereas here `label` is a compact serialized object such
as `"31"` or `"31.210"`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = PartialFlagVariety("31");

julia> dimension(X)
4
```
"""
function PartialFlagVariety(label::AbstractString)
  result = decode_label(String(label))
  mdt = _factors_to_mdt(result.factors)
  PartialFlagVariety(mdt)
end

"""
    zero_locus(label::AbstractString) -> ZeroLocus

Construct a zero locus from a ZeroLocus62 label. The label must contain a
bundle part (i.e. include a `.` separator).

As with [`zero_locus(E)`](@ref), this assumes that the decoded bundle defines a
regular section; the label itself is only a serialization of the ambient
variety and bundle data.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> Z = zero_locus("1.0");

julia> dimension(Z)
0
```
"""
function zero_locus(label::AbstractString)
  result = decode_label(String(label))
  isempty(result.summands) && throw(
    ArgumentError("Label \"$label\" has no bundle part; cannot construct a zero locus.")
  )

  mdt = _factors_to_mdt(result.factors)
  DT = dynkin_type(mdt)
  X = PartialFlagVariety(mdt)

  irr = IrrepLevi[]
  for row in result.summands
    λ = _summand_row_to_weight(row, DT)
    push!(irr, IrrepLevi(mdt, λ))
  end

  E = CompletelyReducibleBundle(X, irr)
  zero_locus(E)
end
