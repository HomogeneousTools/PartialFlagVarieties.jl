export MarkedDynkinType
export marked_nodes, unmarked_nodes, levi_type, levi_rank, central_rank
export is_borel
export central_scaling_factor
export decomposition_matrix, decomposition_matrix_inv
export levi_permutation
export marked_dynkin_diagram

"""
    MarkedDynkinType(DT::Type{<:DynkinType}, marked)

Runtime description of a partial flag variety `G/P`, given by a Dynkin type
`DT` and the tuple of marked simple roots defining the parabolic subgroup `P`.

The marked nodes are stored as runtime data rather than type parameters, while
derived invariants such as the Levi type, decomposition matrix, and dimension
are cached on first use.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> mdt = MarkedDynkinType(TypeA{4}, (2,));

julia> marked_nodes(mdt)
(2,)

julia> unmarked_nodes(mdt)
(1, 3, 4)
```
"""
struct MarkedDynkinType
  dynkin::DataType
  marked::Tuple{Vararg{Int}}

  function MarkedDynkinType(dynkin::DataType, marked::Tuple{Vararg{Int}})
    dynkin <: DynkinType || throw(ArgumentError("Expected a Dynkin type, got $dynkin"))

    R = rank(dynkin)
    for m in marked
      1 <= m <= R || throw(
        ArgumentError(
          "Marked node $m is out of range for $(Lie._type_name(dynkin)) (rank $R)"
        ),
      )
    end

    issorted(marked) || throw(ArgumentError("Marked nodes must be sorted, got $marked"))
    length(unique(marked)) == length(marked) || throw(ArgumentError(
      "Marked nodes must be distinct, got $marked"
    ))

    new(dynkin, marked)
  end
end

function MarkedDynkinType(::Type{DT}, marked::Tuple{Vararg{Int}}) where {DT<:DynkinType}
  marked_int = Tuple(sort(collect(marked)))
  invoke(MarkedDynkinType, Tuple{DataType,Tuple{Vararg{Int}}}, DT, marked_int)
end

function MarkedDynkinType(::Type{DT}, marked::Tuple) where {DT<:DynkinType}
  marked_int = Tuple(sort(Int[m for m in marked]))
  invoke(MarkedDynkinType, Tuple{DataType,Tuple{Vararg{Int}}}, DT, marked_int)
end

MarkedDynkinType(::Type{DT}, marked::Vector{<:Integer}) where {DT<:DynkinType} = MarkedDynkinType(
  DT, Tuple(sort(Int.(marked)))
)

MarkedDynkinType(::Type{DT}, marked::Integer) where {DT<:DynkinType} = MarkedDynkinType(
  DT, (Int(marked),)
)

struct _MarkedDynkinData
  unmarked::Tuple{Vararg{Int}}
  levi::Union{Nothing,DataType}
  levi_permutation::Tuple{Vararg{Int}}
  central_scaling_factor::Int
  decomposition_matrix::Matrix{Rational{Int}}
  decomposition_matrix_inv::Matrix{Rational{Int}}
  dimension::Int
end

const _marked_dynkin_cache = Dict{MarkedDynkinType,_MarkedDynkinData}()

function _compute_marked_dynkin_data(mdt::MarkedDynkinType)
  DT = dynkin_type(mdt)
  marked = marked_nodes(mdt)
  R = rank(DT)

  unmarked = Tuple(i for i in 1:R if !(i in marked))
  if isempty(unmarked)
    levi = nothing
    perm = ()
  else
    C = Lie._cartan_matrix_data(DT)
    idx = collect(unmarked)
    C_sub = C[idx, idx]
    ct = cartan_type(C_sub)
    levi = _cartan_type_to_dynkin_type(ct)
    _, ord = cartan_type_with_ordering(C_sub)
    perm = Tuple(ord)
  end

  Cinv = Lie.cartan_matrix_inverse(DT)
  sf = 1
  for j in marked
    for k in 1:R
      sf = lcm(sf, denominator(Cinv[j, k]))
    end
  end

  M = zeros(Rational{Int}, R, R)
  for i in unmarked
    M[i, i] = 1
  end
  for j in marked
    for k in 1:R
      M[j, k] = Cinv[j, k]
    end
  end

  Minv = inv(M)
  n_pos_G = n_positive_roots(DT)
  n_pos_L = levi === nothing ? 0 : n_positive_roots(levi)

  _MarkedDynkinData(unmarked, levi, perm, sf, M, Minv, n_pos_G - n_pos_L)
end

_mdt_data(mdt::MarkedDynkinType) =
  get!(_marked_dynkin_cache, mdt) do
    _compute_marked_dynkin_data(mdt)
  end

"""Return the ambient Dynkin type of `mdt`."""
dynkin_type(mdt::MarkedDynkinType) = mdt.dynkin

"""Return the marked nodes defining the parabolic subgroup."""
marked_nodes(mdt::MarkedDynkinType) = mdt.marked

"""Return the unmarked simple roots, i.e. the Levi nodes."""
unmarked_nodes(mdt::MarkedDynkinType) = _mdt_data(mdt).unmarked

"""Return the rank of the center of the Levi subgroup."""
central_rank(mdt::MarkedDynkinType) = length(marked_nodes(mdt))

"""Return the Dynkin type of the semisimple Levi factor, or `nothing` for `G/B`."""
levi_type(mdt::MarkedDynkinType) = _mdt_data(mdt).levi

"""Return the rank of the semisimple Levi factor."""
levi_rank(mdt::MarkedDynkinType) = length(unmarked_nodes(mdt))

"""Return `true` exactly for full flag varieties `G/B`."""
is_borel(mdt::MarkedDynkinType) = isempty(unmarked_nodes(mdt))

"""Return the permutation sending natural Levi-node order to canonical Cartan order."""
levi_permutation(mdt::MarkedDynkinType) = _mdt_data(mdt).levi_permutation

"""Return the integer clearing denominators of the marked rows of the decomposition matrix."""
central_scaling_factor(mdt::MarkedDynkinType) = _mdt_data(mdt).central_scaling_factor

"""
    decomposition_matrix(mdt::MarkedDynkinType) -> Matrix{Rational{Int}}

Return the change-of-basis matrix from ambient fundamental weights to
central-plus-Levi coordinates. Unmarked rows are the identity and marked rows
are the corresponding rows of the inverse Cartan matrix.
"""
@inline decomposition_matrix(mdt::MarkedDynkinType) = _mdt_data(mdt).decomposition_matrix

"""Return the inverse of [`decomposition_matrix`](@ref)."""
@inline decomposition_matrix_inv(mdt::MarkedDynkinType) =
  _mdt_data(mdt).decomposition_matrix_inv

Lie.rank(mdt::MarkedDynkinType) = rank(dynkin_type(mdt))

"""Return the dimension of the partial flag variety encoded by `mdt`."""
dimension(mdt::MarkedDynkinType) = _mdt_data(mdt).dimension

function _nonparabolic_height(α_vec::AbstractVector{<:Integer}, marked::Tuple{Vararg{Int}})
  h = 0
  for m in marked
    h += α_vec[m]
  end
  h
end

function positive_nonparabolic_roots(mdt::MarkedDynkinType)
  RS = RootSystem(dynkin_type(mdt))
  roots = collect(positive_roots(RS))
  result = eltype(roots)[]
  for α in roots
    _nonparabolic_height(coefficients(α), marked_nodes(mdt)) > 0 && push!(result, α)
  end
  result
end

function positive_parabolic_roots(mdt::MarkedDynkinType)
  RS = RootSystem(dynkin_type(mdt))
  roots = collect(positive_roots(RS))
  result = eltype(roots)[]
  for α in roots
    _nonparabolic_height(coefficients(α), marked_nodes(mdt)) == 0 && push!(result, α)
  end
  result
end

function tangent_weights(mdt::MarkedDynkinType)
  DT = dynkin_type(mdt)
  RS = RootSystem(DT)
  nonpar_roots = positive_nonparabolic_roots(mdt)
  simple_par = [simple_root(RS, i) for i in unmarked_nodes(mdt)]
  nonpar_set = Set(coefficients(α) for α in nonpar_roots)

  highest = filter(nonpar_roots) do α
    !any(s -> (coefficients(α) + coefficients(s)) in nonpar_set, simple_par)
  end

  [WeightLatticeElem(α) for α in highest]
end

_ambient_type(mdt::MarkedDynkinType) = dynkin_type(mdt)

function Base.show(io::IO, mdt::MarkedDynkinType)
  print(io, "$(Lie._type_name(dynkin_type(mdt))) / P_{$(join(marked_nodes(mdt), ","))}")
end

function marked_dynkin_diagram(mdt::MarkedDynkinType)
  DT = dynkin_type(mdt)
  marked = marked_nodes(mdt)
  diagram = dynkin_diagram(DT)
  lines = split(diagram, '\n')
  marked_set = Set(marked)

  if DT <: SimpleDynkinType
    if DT <: TypeD
      return _marked_diagram_D(DT, marked)
    end
    if DT <: TypeE
      return _marked_diagram_E(DT, marked)
    end

    node_line = lines[1]
    result_chars = collect(node_line)
    node_idx = 0
    for (pos, ch) in enumerate(result_chars)
      if ch == '○'
        node_idx += 1
        if node_idx in marked_set
          result_chars[pos] = '×'
        end
      end
    end
    lines[1] = String(result_chars)
    return join(lines, '\n')
  end

  diagram * "\n(marked: $(join(marked, ", ")))"
end

function _marked_diagram_D(::Type{DT}, marked) where {DT<:TypeD}
  N = rank(DT)
  marked_set = Set(marked)
  node_char(i) = i in marked_set ? '×' : '○'

  prefix = " "^(4 * (N - 2)) * "$(node_char(N)) $N"
  fork = " "^(4 * (N - 2) - 1) * "/"
  if N - 1 >= 2
    main = join([string(node_char(i)) for i in 1:(N - 1)], "───")
    main_labels = join([lpad(string(i), 1) for i in 1:(N - 1)], "   ")
  else
    main = string(node_char(N - 1))
    main_labels = "$(N - 1)"
  end
  prefix * "\n" * fork * "\n" * main * "\n" * main_labels
end

function _marked_diagram_E(::Type{DT}, marked) where {DT<:TypeE}
  N = rank(DT)
  marked_set = Set(marked)
  node_char(i) = i in marked_set ? '×' : '○'

  main_nodes = [1; collect(3:N)]
  main = join([string(node_char(i)) for i in main_nodes], "───")
  main_labels = join([lpad(string(i), 1) for i in main_nodes], "   ")

  indent = 8
  top = " "^indent * "$(node_char(2)) 2"
  branch = " "^indent * "|"
  top * "\n" * branch * "\n" * main * "\n" * main_labels
end

Base.:(==)(a::MarkedDynkinType, b::MarkedDynkinType) =
  dynkin_type(a) == dynkin_type(b) && marked_nodes(a) == marked_nodes(b)

Base.hash(mdt::MarkedDynkinType, h::UInt) = hash((dynkin_type(mdt), marked_nodes(mdt)), h)
