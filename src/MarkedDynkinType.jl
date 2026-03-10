# ═══════════════════════════════════════════════════════════════════════════════
#  MarkedDynkinType — a Dynkin type with crossed-out nodes (type-level encoding)
#
#  A marked Dynkin diagram encodes both the partial flag variety G/P_I
#  (where I is the set of crossed-out / nonparabolic / removed nodes)
#  and the Levi factor of P_I: the unmarked nodes determine the semisimple
#  part of the Levi subgroup, while the marked nodes correspond to the
#  center of the Levi (whose rank equals the Picard rank of G/P_I).
#
#  Both the Dynkin type and the set of marked nodes are encoded as type
#  parameters, enabling heavy compile-time specialization following the
#  same pattern as Lie.jl.
# ═══════════════════════════════════════════════════════════════════════════════

export MarkedDynkinType
export marked_nodes, unmarked_nodes, levi_type, levi_rank, central_rank
export central_scaling_factor
export decomposition_matrix, decomposition_matrix_inv
export levi_permutation
export marked_dynkin_diagram

# Names from Lie, StaticArrays, LinearAlgebra are available via the parent module.

# ═══════════════════════════════════════════════════════════════════════════════
#  Core type
# ═══════════════════════════════════════════════════════════════════════════════

"""
    MarkedDynkinType{DT, Marked}

A Dynkin type `DT` with a set of **marked** (crossed-out / nonparabolic) nodes
encoded as a compile-time sorted tuple `Marked`.

The marked nodes are the nodes removed from the Dynkin diagram to define
the parabolic subgroup ``P_I``, so that the partial flag variety is ``G/P_I``.
The unmarked nodes determine the Levi factor of ``P_I``: they give the
semisimple part of the Levi subgroup, while the marked nodes correspond
to its center.

# Type parameters
- `DT <: DynkinType`: the ambient Dynkin type (e.g., `TypeA{4}`)
- `Marked`: a sorted `Tuple` of `Int` indices of the crossed-out nodes

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> MDT = MarkedDynkinType{TypeA{4}, (2,)}
MarkedDynkinType{A4, (2,)}

julia> marked_nodes(MDT)
(2,)

julia> unmarked_nodes(MDT)
(1, 3, 4)

julia> levi_rank(MDT)
3

julia> central_rank(MDT)
1
```
"""
struct MarkedDynkinType{DT<:DynkinType,Marked} end

function MarkedDynkinType{DT,Marked}(::Nothing) where {DT<:DynkinType,Marked}
  R = rank(DT)
  for m in Marked
    1 <= m <= R || throw(ArgumentError(
      "Marked node $m is out of range for $(Lie._type_name(DT)) (rank $R)"
    ))
  end
  issorted(Marked) || throw(ArgumentError(
    "Marked nodes must be sorted, got $Marked"
  ))
  length(unique(Marked)) == length(Marked) || throw(ArgumentError(
    "Marked nodes must be distinct, got $Marked"
  ))
  return MarkedDynkinType{DT,Marked}()
end

"""
    MarkedDynkinType(::Type{DT}, marked::NTuple{K,Int}) where {DT, K}

Construct a `MarkedDynkinType` from a Dynkin type and a tuple of marked node indices.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> MarkedDynkinType(TypeA{3}, (2,))
A3 / P_{2}

julia> MarkedDynkinType(TypeB{4}, (1, 3))
B4 / P_{1,3}
```
"""
function MarkedDynkinType(::Type{DT}, marked::NTuple{K,Int}) where {DT<:DynkinType,K}
  sorted = Tuple(sort(collect(marked)))
  return MarkedDynkinType{DT,sorted}(nothing)
end

function MarkedDynkinType(::Type{DT}, marked::Vector{<:Integer}) where {DT<:DynkinType}
  return MarkedDynkinType(DT, Tuple(sort(Int.(marked))))
end

function MarkedDynkinType(::Type{DT}, marked::Integer) where {DT<:DynkinType}
  return MarkedDynkinType(DT, (Int(marked),))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Basic accessors (all compile-time)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    marked_nodes(::Type{MarkedDynkinType{DT,Marked}}) -> Tuple

Return the indices of the crossed-out (nonparabolic) nodes.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> marked_nodes(MarkedDynkinType{TypeA{4}, (2,)})
(2,)

julia> marked_nodes(MarkedDynkinType{TypeD{5}, (1, 5)})
(1, 5)
```
"""
marked_nodes(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked} = Marked

"""
    unmarked_nodes(::Type{MarkedDynkinType{DT,Marked}}) -> Tuple

Return the indices of the parabolic (Levi / kept) nodes.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> unmarked_nodes(MarkedDynkinType{TypeA{4}, (2,)})
(1, 3, 4)
```
"""
@generated function unmarked_nodes(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  result = Tuple(i for i in 1:R if !(i in Marked))
  return :($result)
end

"""
    central_rank(::Type{MDT}) -> Int

Return the number of marked (crossed-out) nodes, which equals the Picard rank
of the partial flag variety ``G/P``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> central_rank(MarkedDynkinType{TypeA{4}, (2,)})
1

julia> central_rank(MarkedDynkinType{TypeA{3}, (1, 3)})
2
```
"""
central_rank(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked} = length(Marked)

Lie.rank(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked} = rank(DT)

"""
    central_scaling_factor(::Type{MDT}) -> Int

Return the scaling factor for the central part of Levi representations
associated to `MDT`.  This is the LCM of all denominators of the
inverse Cartan matrix entries at the marked rows:

``\\mathrm{lcm}\\{\\mathrm{denom}(C^{-1}[j, k]) : j \\in \\mathrm{Marked},\\; 1 \\le k \\le R\\}``

Central characters are stored internally as integers multiplied by this
factor, eliminating all `Rational{Int}` arithmetic from hot paths.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> central_scaling_factor(MarkedDynkinType{TypeA{3}, (2,)})
2

julia> central_scaling_factor(MarkedDynkinType{TypeA{4}, (2,)})
5
```
"""
@generated function central_scaling_factor(
  ::Type{MDT},
) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  R = rank(DT)
  Cinv = Lie.cartan_matrix_inverse(DT)
  sf = 1
  for j in Marked
    for k in 1:R
      sf = lcm(sf, denominator(Cinv[j, k]))
    end
  end
  return :($sf)
end

central_scaling_factor(mdt::MarkedDynkinType) = central_scaling_factor(typeof(mdt))

# ═══════════════════════════════════════════════════════════════════════════════
#  Levi type computation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    levi_type(::Type{MarkedDynkinType{DT,Marked}}) -> Type{<:DynkinType}

Return the Dynkin type of the semisimple part of the Levi subgroup.
This is computed from the sub-Cartan matrix on the unmarked nodes,
using `cartan_type` to identify the resulting Cartan matrix.

Returns `nothing` for the full flag variety (all nodes marked), since the
Levi is then a torus with trivial semisimple part.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> levi_type(MarkedDynkinType{TypeA{4}, (2,)})
ProductDynkinType{Tuple{TypeA{1}, TypeA{2}}}

julia> levi_type(MarkedDynkinType{TypeA{4}, (3,)})
ProductDynkinType{Tuple{TypeA{2}, TypeA{1}}}

julia> levi_type(MarkedDynkinType{TypeD{5}, (5,)})
TypeA{4}

julia> levi_type(MarkedDynkinType{TypeA{3}, (1, 2, 3)}) === nothing
true
```
"""
@generated function levi_type(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  unmarked = [i for i in 1:R if !(i in Marked)]
  n = length(unmarked)
  if n == 0
    return :(nothing)
  end
  C = Lie._cartan_matrix_data(DT)
  C_sub = C[unmarked, unmarked]
  ct = cartan_type(C_sub)
  lt = _cartan_type_to_dynkin_type(ct)
  return :($lt)
end

# Instance-dispatch overload: allows `levi_type(mdt)` where `mdt` is a
# MarkedDynkinType value (as returned by `marked_dynkin_type(rep)`).
levi_type(mdt::MarkedDynkinType) = levi_type(typeof(mdt))

"""
    levi_permutation(::Type{MDT}) -> Tuple{Int,...}

Return the permutation that converts from natural (sub-diagram position) ordering
of the unmarked nodes to the canonical Dynkin ordering of `levi_type(MDT)`.

For `i` in `1:levi_rank(MDT)`, the `i`-th fundamental weight of `levi_type` corresponds
to the `perm[i]`-th unmarked node (in 1-based position order of unmarked nodes).

This permutation is computed via `cartan_type_with_ordering` of the sub-Cartan matrix.
"""
@generated function levi_permutation(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  unmarked = [i for i in 1:R if !(i in Marked)]
  n = length(unmarked)
  if n == 0
    return :(())
  end
  C = Lie._cartan_matrix_data(DT)
  C_sub = C[unmarked, unmarked]
  _, ord = cartan_type_with_ordering(C_sub)
  # ord[i] = which sub-node corresponds to canonical LT node i
  # i.e., ss_coords_LT[i] <- ss_coords_nat[ord[i]]
  perm_tuple = Tuple(ord)
  return :($perm_tuple)
end

"""
    levi_rank(::Type{MDT}) -> Int

Return the rank of the semisimple part of the Levi subgroup.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> levi_rank(MarkedDynkinType{TypeA{4}, (2,)})
3

julia> levi_rank(MarkedDynkinType{TypeE{6}, (1, 6)})
4
```
"""
@generated function levi_rank(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  n = R - length(Marked)
  return :($n)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Decomposition matrix (change of basis)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    decomposition_matrix(::Type{MarkedDynkinType{DT,Marked}}) -> SMatrix

Return the change-of-basis matrix ``M`` from fundamental weight coordinates
to the (central + semisimple) coordinate system of the Levi decomposition.

Given a weight ``\\lambda = \\sum_i \\lambda_i \\omega_i``, the new coordinates
``\\lambda' = M \\lambda`` split as:
- ``\\lambda'_j`` for ``j \\in`` marked (nonparabolic) nodes: **central part**
  (character of the center of the Levi)
- ``\\lambda'_i`` for ``i \\in`` unmarked (parabolic) nodes: **semisimple part**
  (weight of the semisimple part of the Levi)

The matrix is constructed as:
- Row ``i`` for ``i`` unmarked (parabolic): standard basis vector ``e_i``
- Row ``j`` for ``j`` marked (nonparabolic): ``C^{-1}[j, :]``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> M = decomposition_matrix(MarkedDynkinType{TypeA{3}, (2,)});

julia> M[1, 1] == 1  # parabolic row
true

julia> M[3, 3] == 1  # parabolic row
true
```
"""
@generated function decomposition_matrix(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  Cinv = Lie.cartan_matrix_inverse(DT)

  unmarked = [i for i in 1:R if !(i in Marked)]

  M = zeros(Rational{Int}, R, R)
  # Parabolic (unmarked) rows: identity
  for i in unmarked
    M[i, i] = 1
  end
  # Nonparabolic (marked) rows: C^{-1}[j, :]
  for j in Marked
    for k in 1:R
      M[j, k] = Cinv[j, k]
    end
  end

  entries = Tuple(M[i, j] for j in 1:R for i in 1:R)
  return :(SMatrix{$R, $R, Rational{Int}, $(R * R)}($entries))
end

"""
    decomposition_matrix_inv(::Type{MarkedDynkinType{DT,Marked}}) -> SMatrix

Return the inverse of [`decomposition_matrix`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties, StaticArrays

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> M = decomposition_matrix(MDT);

julia> Minv = decomposition_matrix_inv(MDT);

julia> using StaticArrays: SMatrix

julia> using LinearAlgebra: I

julia> M * Minv ≈ SMatrix{3,3}(1.0I)
true
```
"""
@generated function decomposition_matrix_inv(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  Cinv = Lie.cartan_matrix_inverse(DT)

  unmarked = [i for i in 1:R if !(i in Marked)]

  M = zeros(Rational{Int}, R, R)
  for i in unmarked
    M[i, i] = 1
  end
  for j in Marked
    for k in 1:R
      M[j, k] = Cinv[j, k]
    end
  end

  Minv = inv(M)
  entries = Tuple(Minv[i, j] for j in 1:R for i in 1:R)
  return :(SMatrix{$R, $R, Rational{Int}, $(R * R)}($entries))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Positive roots decomposition (parabolic vs non-parabolic)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _nonparabolic_height(α::RootSpaceElem, marked::Tuple) -> Int

Sum of the coefficients of `α` at the marked (nonparabolic) node positions.
"""
function _nonparabolic_height(α_vec::SVector{R,Int}, marked::NTuple{K,Int}) where {R,K}
  h = 0
  for m in marked
    h += α_vec[m]
  end
  return h
end

"""
    positive_nonparabolic_roots(::Type{MDT}) -> Vector{RootSpaceElem}

Return the positive roots that have nonzero coefficient at some marked node.
These are the roots of the unipotent radical of `P`.
"""
function positive_nonparabolic_roots(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  RS = RootSystem(DT)
  R = rank(DT)
  result = RootSpaceElem{DT,R}[]
  for α in positive_roots(RS)
    _nonparabolic_height(coefficients(α), Marked) > 0 && push!(result, α)
  end
  return result
end

"""
    positive_parabolic_roots(::Type{MDT}) -> Vector{RootSpaceElem}

Return the positive roots that have zero coefficient at all marked nodes.
These are the positive roots of the Levi subgroup.
"""
function positive_parabolic_roots(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  RS = RootSystem(DT)
  R = rank(DT)
  result = RootSpaceElem{DT,R}[]
  for α in positive_roots(RS)
    _nonparabolic_height(coefficients(α), Marked) == 0 && push!(result, α)
  end
  return result
end

"""
    tangent_weights(::Type{MDT}) -> Vector{WeightLatticeElem}

Return the highest weights of the tangent bundle ``T_{G/P}`` viewed as a
completely reducible homogeneous bundle (its semisimplification).

These are the positive nonparabolic roots ``\\alpha`` such that
``\\alpha + \\beta`` is **not** a positive nonparabolic root for any simple
parabolic root ``\\beta`` (i.e., ``\\alpha`` is maximal in the poset of
nonparabolic roots under the Levi root ordering).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> length(tangent_weights(MarkedDynkinType{TypeA{3}, (2,)}))  # Gr(2,4): T = S* ⊗ Q (irreducible)
1

julia> length(tangent_weights(MarkedDynkinType{TypeA{4}, (1,)}))  # ℙ⁴
1
```
"""
function tangent_weights(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  R = rank(DT)
  RS = RootSystem(DT)
  Marked = marked_nodes(MDT)
  unmarked = unmarked_nodes(MDT)

  nonpar_roots = positive_nonparabolic_roots(MDT)

  # Simple parabolic roots (simple roots at unmarked positions)
  simple_par = [simple_root(RS, i) for i in unmarked]

  # Set of nonparabolic root vectors for fast lookup
  nonpar_set = Set(coefficients(α) for α in nonpar_roots)

  # Filter: keep roots where adding any simple parabolic root
  # gives something NOT in the nonparabolic root set
  highest = filter(nonpar_roots) do α
    !any(s -> (coefficients(α) + coefficients(s)) in nonpar_set, simple_par)
  end

  # Convert to weight lattice elements
  return [WeightLatticeElem(α) for α in highest]
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Helper: extract ambient DynkinType from MDT
# ═══════════════════════════════════════════════════════════════════════════════

_ambient_type(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked} = DT

# ═══════════════════════════════════════════════════════════════════════════════
#  Pretty-printing
# ═══════════════════════════════════════════════════════════════════════════════

function Base.show(io::IO, ::MarkedDynkinType{DT,Marked}) where {DT,Marked}
  print(io, "$(Lie._type_name(DT)) / P_{$(join(Marked, ","))}")
end

function Base.show(io::IO, ::Type{MarkedDynkinType{DT,Marked}}) where {DT<:DynkinType,Marked}
  print(io, "MarkedDynkinType{$(Lie._type_name(DT)), $Marked}")
end

"""
    marked_dynkin_diagram(::Type{MDT}) -> String

Render the Dynkin diagram with marked (crossed-out) nodes shown as `×`
and unmarked (parabolic) nodes shown as `○`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> println(marked_dynkin_diagram(MarkedDynkinType{TypeA{4}, (2,)}))
○───×───○───○
1   2   3   4

julia> println(marked_dynkin_diagram(MarkedDynkinType{TypeB{3}, (1,)}))
×───○═>═○
1   2   3
```
"""
function marked_dynkin_diagram(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  diagram = dynkin_diagram(DT)
  lines = split(diagram, '\n')

  marked_set = Set(Marked)
  R = rank(DT)

  if DT <: SimpleDynkinType
    if DT <: TypeD
      return _marked_diagram_D(DT, Marked)
    end
    if DT <: TypeE
      return _marked_diagram_E(DT, Marked)
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

  return diagram * "\n(marked: $(join(Marked, ", ")))"
end

function _marked_diagram_D(::Type{DT}, Marked) where {DT<:TypeD}
  N = rank(DT)
  marked_set = Set(Marked)
  node_char(i) = i in marked_set ? '×' : '○'

  prefix = " "^(4 * (N - 2)) * "$(node_char(N)) $N"
  fork = " "^(4 * (N - 2) - 1) * "/"
  if N - 1 >= 2
    main = join([string(node_char(i)) for i in 1:(N - 1)], "───")
    main_labels = join([lpad(string(i), 1) for i in 1:(N - 1)], "   ")
  else
    main = string(node_char(N - 1))
    main_labels = "$(N-1)"
  end
  return prefix * "\n" * fork * "\n" * main * "\n" * main_labels
end

function _marked_diagram_E(::Type{DT}, Marked) where {DT<:TypeE}
  N = rank(DT)
  marked_set = Set(Marked)
  node_char(i) = i in marked_set ? '×' : '○'

  n_main = N - 1
  main_nodes = [1; collect(3:N)]
  main = join([string(node_char(i)) for i in main_nodes], "───")
  main_labels = join([lpad(string(i), 1) for i in main_nodes], "   ")

  indent = 8
  top = " "^indent * "$(node_char(2)) 2"
  branch = " "^indent * "|"
  return top * "\n" * branch * "\n" * main * "\n" * main_labels
end
