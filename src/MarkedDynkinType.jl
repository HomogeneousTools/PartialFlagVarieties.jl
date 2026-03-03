# ═══════════════════════════════════════════════════════════════════════════════
#  MarkedDynkinType — a Dynkin type with crossed-out nodes (type-level encoding)
#
#  A marked Dynkin diagram encodes a partial flag variety G/P_I where I is
#  the set of crossed-out (nonparabolic / removed) nodes. The remaining
#  (unmarked) nodes determine the Levi subgroup.
#
#  Both the Dynkin type and the set of marked nodes are encoded as type
#  parameters, enabling heavy compile-time specialization following the
#  same pattern as Lie.jl.
# ═══════════════════════════════════════════════════════════════════════════════

export MarkedDynkinType
export marked_nodes, unmarked_nodes, levi_type, levi_rank, central_rank
export special_matrix, special_matrix_inv
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
The unmarked nodes determine the Levi factor of ``P_I``.

# Type parameters
- `DT <: DynkinType`: the ambient Dynkin type (e.g., `TypeA{4}`)
- `Marked`: a sorted `Tuple` of `Int` indices of the crossed-out nodes

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> MDT = MarkedDynkinType{TypeA{4}, (2,)}
MarkedDynkinType{TypeA{4}, (2,)}

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
MarkedDynkinType{TypeA{3}, (2,)}

julia> MarkedDynkinType(TypeB{4}, (1, 3))
MarkedDynkinType{TypeB{4}, (1, 3)}
```
"""
function MarkedDynkinType(::Type{DT}, marked::NTuple{K,Int}) where {DT<:DynkinType,K}
  sorted = Tuple(sort(collect(marked)))
  return MarkedDynkinType{DT,sorted}(nothing)
end

function MarkedDynkinType(::Type{DT}, marked::Vector{Int}) where {DT<:DynkinType}
  return MarkedDynkinType(DT, Tuple(sort(marked)))
end

function MarkedDynkinType(::Type{DT}, marked::Int) where {DT<:DynkinType}
  return MarkedDynkinType(DT, (marked,))
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

# ═══════════════════════════════════════════════════════════════════════════════
#  Levi type computation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _levi_dynkin_type(C_sub::Matrix{Int}) -> Type{<:DynkinType}

Given a sub-Cartan matrix (for the unmarked nodes), determine the
`DynkinType` of the semisimple part of the Levi subgroup by decomposing
into connected components and identifying each component.
"""
function _levi_dynkin_type(C_sub::AbstractMatrix{Int}, n::Int)
  if n == 0
    # Full flag variety: Levi is a torus, semisimple part is trivial
    # Return TypeA{1} as a dummy (rank 0 isn't representable);
    # we handle this specially
    return nothing
  end

  # Find connected components via adjacency in the Cartan matrix
  adj = [Int[] for _ in 1:n]
  for i in 1:n, j in (i + 1):n
    if C_sub[i, j] != 0 || C_sub[j, i] != 0
      push!(adj[i], j)
      push!(adj[j], i)
    end
  end

  visited = falses(n)
  components = Vector{Vector{Int}}()
  for start in 1:n
    visited[start] && continue
    comp = Int[]
    stack = [start]
    visited[start] = true
    while !isempty(stack)
      u = pop!(stack)
      push!(comp, u)
      for v in adj[u]
        if !visited[v]
          visited[v] = true
          push!(stack, v)
        end
      end
    end
    sort!(comp)
    push!(components, comp)
  end

  # Identify each connected component
  component_types = []
  for comp in components
    C_comp = C_sub[comp, comp]
    r = length(comp)
    dt = _identify_cartan_matrix(C_comp, r)
    push!(component_types, dt)
  end

  if length(component_types) == 1
    return component_types[1]
  else
    # Build ProductDynkinType
    return _make_product_type(component_types)
  end
end

"""
    _identify_cartan_matrix(C::AbstractMatrix{Int}, r::Int) -> Type{<:SimpleDynkinType}

Identify a connected Cartan matrix of rank `r`, handling arbitrary
row/column orderings (permuted sub-Cartan matrices).

Uses graph structure: builds the adjacency graph, then classifies by
degree sequence and branch lengths.
"""
function _identify_cartan_matrix(C::AbstractMatrix{Int}, r::Int)
  # A_1
  r == 1 && return TypeA{1}

  # Build adjacency graph and detect edge types
  adj = [Int[] for _ in 1:r]
  has_double = false
  has_triple = false
  double_edge = (0, 0)  # (i, j) where C[i,j] = -2 (i is the "short" side)
  for i in 1:r, j in (i+1):r
    C[i, j] == 0 && C[j, i] == 0 && continue
    push!(adj[i], j)
    push!(adj[j], i)
    if C[i, j] == -2
      has_double = true
      double_edge = (i, j)
    elseif C[j, i] == -2
      has_double = true
      double_edge = (j, i)
    end
    if C[i, j] == -3 || C[j, i] == -3
      has_triple = true
    end
  end

  degrees = [length(adj[i]) for i in 1:r]
  max_deg = maximum(degrees)

  # ── Simply-laced ──────────────────────────────────────────────────────
  if !has_double && !has_triple
    if max_deg <= 2
      # Path graph → A_r
      return TypeA{r}
    end

    # Has a node of degree 3 → D_r or E_r
    branch_node = findfirst(==(3), degrees)
    branch_node === nothing && error(
      "Could not identify simply-laced Cartan matrix of rank $r: $C"
    )

    # Compute branch lengths (distance from branch node to leaf in each direction)
    branch_lengths = Int[]
    for neighbor in adj[branch_node]
      len = 1
      prev = branch_node
      curr = neighbor
      while length(adj[curr]) == 2
        next = (adj[curr][1] == prev) ? adj[curr][2] : adj[curr][1]
        prev = curr
        curr = next
        len += 1
      end
      push!(branch_lengths, len)
    end
    sort!(branch_lengths)

    # D_r: branch lengths [1, 1, r-3]
    if branch_lengths == [1, 1, r - 3]
      return TypeD{r}
    end
    # E_6: [1, 2, 2], E_7: [1, 2, 3], E_8: [1, 2, 4]
    if branch_lengths == [1, 2, 2] && r == 6
      return TypeE{6}
    elseif branch_lengths == [1, 2, 3] && r == 7
      return TypeE{7}
    elseif branch_lengths == [1, 2, 4] && r == 8
      return TypeE{8}
    end

    error("Could not identify simply-laced Cartan matrix of rank $r " *
          "with branch lengths $branch_lengths: $C")
  end

  # ── Non-simply-laced ──────────────────────────────────────────────────

  # G2: rank 2, triple bond
  if has_triple && r == 2
    return TypeG2
  end

  # B2/C2: rank 2, double bond
  if has_double && r == 2
    # B2 and C2 are isomorphic; Lie.jl uses B2
    return TypeB{2}
  end

  # F4: rank 4, double bond in the middle of a path
  if has_double && r == 4 && max_deg <= 2
    (di, dj) = double_edge
    # Both sides of the double bond should be interior (degree 2)
    if degrees[di] == 2 && degrees[dj] == 2
      return TypeF4
    end
  end

  # B_r or C_r: path with double bond at one end
  if has_double && max_deg <= 2
    (di, dj) = double_edge  # C[di, dj] = -2
    # di is the "short root" side
    # B_r: short root is a leaf → di is a leaf
    if degrees[di] == 1
      return TypeB{r}
    end
    # C_r: short root is not a leaf → di is interior, dj is leaf
    if degrees[dj] == 1
      return TypeC{r}
    end
  end

  error("Could not identify Cartan matrix of rank $r: $C")
end

function _make_product_type(types::Vector)
  if length(types) == 1
    return types[1]
  end
  # Build the ProductDynkinType at the type level
  # We need to construct: ProductDynkinType{Tuple{T1, T2, ...}}
  return ProductDynkinType{Tuple{types...}}
end

"""
    levi_type(::Type{MarkedDynkinType{DT,Marked}}) -> Type{<:DynkinType}

Return the Dynkin type of the semisimple part of the Levi subgroup.
This is computed from the sub-Cartan matrix on the unmarked nodes.

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
  um = [i for i in 1:R if !(i in Marked)]
  n = length(um)
  if n == 0
    return :(nothing)
  end
  C = Lie._cartan_matrix_data(DT)
  C_sub = C[um, um]
  lt = _levi_dynkin_type(C_sub, n)
  return :($lt)
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
#  Special matrix (change of basis)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    special_matrix(::Type{MarkedDynkinType{DT,Marked}}) -> SMatrix

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

julia> M = special_matrix(MarkedDynkinType{TypeA{3}, (2,)});

julia> M[1, 1] == 1  # parabolic row
true

julia> M[3, 3] == 1  # parabolic row
true
```
"""
@generated function special_matrix(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  C = Lie._cartan_matrix_data(DT)
  Crat = Rational{Int}.(C)
  Cinv = inv(Crat)

  um = [i for i in 1:R if !(i in Marked)]

  M = zeros(Rational{Int}, R, R)
  # Parabolic (unmarked) rows: identity
  for i in um
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
    special_matrix_inv(::Type{MarkedDynkinType{DT,Marked}}) -> SMatrix

Return the inverse of [`special_matrix`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties, StaticArrays

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> M = special_matrix(MDT);

julia> Minv = special_matrix_inv(MDT);

julia> M * Minv ≈ SMatrix{3,3}(1.0I)
true
```
"""
@generated function special_matrix_inv(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  C = Lie._cartan_matrix_data(DT)
  Crat = Rational{Int}.(C)
  Cinv = inv(Crat)

  um = [i for i in 1:R if !(i in Marked)]

  M = zeros(Rational{Int}, R, R)
  for i in um
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
#  Dimension and topological invariants
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dimension(::Type{MDT}) -> Int
    dimension(::Type{MDT}) -> Int

Return the dimension of the partial flag variety ``G/P``.

This equals the number of positive roots of ``G`` that are not roots of the Levi:
``\\dim(G/P) = |\\Phi^+_G| - |\\Phi^+_L|``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> dimension(MarkedDynkinType{TypeA{3}, (2,)})  # Gr(2,4)
4

julia> dimension(MarkedDynkinType{TypeA{4}, (1,)})  # ℙ⁴
4

julia> dimension(MarkedDynkinType{TypeE{6}, (1,)})  # Cayley plane
16
```
"""
@generated function dimension(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  um = [i for i in 1:R if !(i in Marked)]
  n_pos_G = n_positive_roots(DT)
  if isempty(um)
    n_pos_L = 0
  else
    C = Lie._cartan_matrix_data(DT)
    C_sub = C[um, um]
    lt = _levi_dynkin_type(C_sub, length(um))
    n_pos_L = lt === nothing ? 0 : n_positive_roots(lt)
  end
  d = n_pos_G - n_pos_L
  return :($d)
end

"""
    picard_rank(::Type{MDT}) -> Int

Return the Picard rank of ``G/P``, which equals the number of marked nodes.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> picard_rank(MarkedDynkinType{TypeA{3}, (2,)})
1

julia> picard_rank(MarkedDynkinType{TypeA{3}, (1, 3)})
2
```
"""
picard_rank(::Type{MDT}) where {MDT<:MarkedDynkinType} = central_rank(MDT)

"""
    euler_characteristic(::Type{MDT}) -> BigInt

Return the Euler characteristic ``\\chi(G/P) = |W_G| / |W_L|``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> euler_characteristic(MarkedDynkinType{TypeA{3}, (2,)})  # Gr(2,4) = 6
6

julia> euler_characteristic(MarkedDynkinType{TypeA{4}, (1,)})  # ℙ⁴ = 5
5

julia> euler_characteristic(MarkedDynkinType{TypeE{6}, (1,)})  # OP²
27
```
"""
@generated function euler_characteristic(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  um = [i for i in 1:R if !(i in Marked)]
  wG = weyl_order(DT)
  if isempty(um)
    wL = BigInt(1)
  else
    C = Lie._cartan_matrix_data(DT)
    C_sub = C[um, um]
    lt = _levi_dynkin_type(C_sub, length(um))
    wL = lt === nothing ? BigInt(1) : weyl_order(lt)
  end
  chi = wG ÷ wL
  return :($chi)
end

"""
    betti_numbers(::Type{MDT}) -> Vector{BigInt}

Compute the Betti numbers of ``G/P`` using the Poincaré polynomial formula:

``P(t) = \\prod_i \\frac{1 - t^{d_i^G}}{1 - t} \\bigg/ \\prod_j \\frac{1 - t^{d_j^L}}{1 - t}``

where ``d_i^G`` and ``d_j^L`` are the degrees of the fundamental invariants of
the Weyl groups of ``G`` and ``L``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> betti_numbers(MarkedDynkinType{TypeA{2}, (1,)})  # ℙ²
3-element Vector{BigInt}:
 1
 1
 1

julia> betti_numbers(MarkedDynkinType{TypeA{3}, (2,)})  # Gr(2,4)
5-element Vector{BigInt}:
 1
 1
 2
 1
 1
```
"""
@generated function betti_numbers(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  um = [i for i in 1:R if !(i in Marked)]

  degs_G = collect(degrees_fundamental_invariants(DT))

  if isempty(um)
    degs_L = Int[]
  else
    C = Lie._cartan_matrix_data(DT)
    C_sub = C[um, um]
    lt = _levi_dynkin_type(C_sub, length(um))
    degs_L = lt === nothing ? Int[] : collect(degrees_fundamental_invariants(lt))
  end

  # Compute Poincaré polynomial as a coefficient vector
  # P(t) = ∏(1 + t + ... + t^{d-1}) for each d in degs_G
  #       / ∏(1 + t + ... + t^{d-1}) for each d in degs_L
  # Compute numerator polynomial
  function expand_factor(d)
    return ones(BigInt, d)  # [1, 1, ..., 1] of length d = 1 + t + ... + t^{d-1}
  end

  function poly_mul(a, b)
    la, lb = length(a), length(b)
    result = zeros(BigInt, la + lb - 1)
    for i in 1:la, j in 1:lb
      result[i + j - 1] += a[i] * b[j]
    end
    return result
  end

  function poly_div(a, b)
    # Exact polynomial division
    la, lb = length(a), length(b)
    a = copy(a)
    result = zeros(BigInt, la - lb + 1)
    for i in (la - lb + 1):-1:1
      c = a[i + lb - 1] ÷ b[lb]
      result[i] = c
      for j in 1:lb
        a[i + j - 1] -= c * b[j]
      end
    end
    return result
  end

  num = BigInt[1]
  for d in degs_G
    num = poly_mul(num, expand_factor(d))
  end

  den = BigInt[1]
  for d in degs_L
    den = poly_mul(den, expand_factor(d))
  end

  result = poly_div(num, den)

  return :($(Vector{BigInt}(result)))
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
  um = unmarked_nodes(MDT)

  nonpar_roots = positive_nonparabolic_roots(MDT)

  # Simple parabolic roots (simple roots at unmarked positions)
  simple_par = [simple_root(RS, i) for i in um]

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
#  Classification predicates
# ═══════════════════════════════════════════════════════════════════════════════

"""
    is_generalized_grassmannian(::Type{MDT}) -> Bool

A partial flag variety is a generalized Grassmannian iff exactly one node is marked.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_generalized_grassmannian(MarkedDynkinType{TypeA{3}, (2,)})
true

julia> is_generalized_grassmannian(MarkedDynkinType{TypeA{3}, (1, 3)})
false
```
"""
is_generalized_grassmannian(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked} =
  length(Marked) == 1

"""
    is_full_flag(::Type{MDT}) -> Bool

Check if all nodes are marked (i.e., ``G/B``).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_full_flag(MarkedDynkinType{TypeA{2}, (1, 2)})
true

julia> is_full_flag(MarkedDynkinType{TypeA{3}, (2,)})
false
```
"""
is_full_flag(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked} =
  length(Marked) == rank(DT)

"""
    is_cominuscule(::Type{MDT}) -> Bool

Check whether the partial flag variety is cominuscule (each irreducible factor
is a cominuscule generalized Grassmannian).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_cominuscule(MarkedDynkinType{TypeA{3}, (2,)})  # Gr(2,4)
true

julia> is_cominuscule(MarkedDynkinType{TypeB{3}, (2,)})
false

julia> is_cominuscule(MarkedDynkinType{TypeB{3}, (1,)})
true
```
"""
function is_cominuscule(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)

  # Must be a generalized Grassmannian for simple types
  # For product types we'd need to check factors; for simple types:
  DT <: SimpleDynkinType || return false  # TODO: handle products
  length(Marked) != 1 && return false

  m = Marked[1]
  R = rank(DT)

  DT <: TypeA && return true
  DT <: TypeB && return m == 1
  DT <: TypeC && return m == R
  DT <: TypeD && return (m == 1 || m == R - 1 || m == R)
  DT <: TypeE{6} && return (m == 1 || m == 6)
  DT <: TypeE{7} && return m == 7
  return false
end

"""
    is_minuscule(::Type{MDT}) -> Bool

Check whether the partial flag variety is minuscule.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_minuscule(MarkedDynkinType{TypeA{3}, (2,)})
true

julia> is_minuscule(MarkedDynkinType{TypeB{3}, (3,)})
true

julia> is_minuscule(MarkedDynkinType{TypeB{3}, (1,)})
false
```
"""
function is_minuscule(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)

  DT <: SimpleDynkinType || return false
  length(Marked) != 1 && return false

  m = Marked[1]
  R = rank(DT)

  DT <: TypeA && return true
  DT <: TypeB && return m == R
  DT <: TypeC && return m == 1
  DT <: TypeD && return (m == 1 || m == R - 1 || m == R)
  DT <: TypeE{6} && return (m == 1 || m == 6)
  DT <: TypeE{7} && return m == 7
  return false
end

"""
    is_adjoint(::Type{MDT}) -> Bool

Check whether the partial flag variety is the adjoint variety.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_adjoint(MarkedDynkinType{TypeB{3}, (2,)})
true

julia> is_adjoint(MarkedDynkinType{TypeG2, (2,)})
true

julia> is_adjoint(MarkedDynkinType{TypeA{3}, (1, 3)})
true
```
"""
function is_adjoint(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  R = rank(DT)

  DT <: SimpleDynkinType || return false

  # Type A special case: adjoint has nodes 1 and r marked
  if DT <: TypeA
    return length(Marked) == 2 && Marked == (1, R)
  end

  length(Marked) != 1 && return false
  m = Marked[1]

  DT <: TypeB && return m == 2
  DT <: TypeC && return m == 1
  DT <: TypeD && R >= 4 && return m == 2
  DT <: TypeE{6} && return m == 2
  DT <: TypeE{7} && return m == 1
  DT <: TypeE{8} && return m == 8
  DT <: TypeF4 && return m == 1
  DT <: TypeG2 && return m == 2
  return false
end

"""
    is_coadjoint(::Type{MDT}) -> Bool

Check whether the partial flag variety is the coadjoint variety.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_coadjoint(MarkedDynkinType{TypeB{3}, (1,)})
true

julia> is_coadjoint(MarkedDynkinType{TypeG2, (1,)})
true
```
"""
function is_coadjoint(::Type{MDT}) where {MDT<:MarkedDynkinType}
  DT = _ambient_type(MDT)
  Marked = marked_nodes(MDT)
  R = rank(DT)

  DT <: SimpleDynkinType || return false

  # Type A: same as adjoint (self-dual)
  if DT <: TypeA
    return length(Marked) == 2 && Marked == (1, R)
  end

  length(Marked) != 1 && return false
  m = Marked[1]

  DT <: TypeB && return m == 1
  DT <: TypeC && return m == 2
  DT <: TypeD && R >= 4 && return m == 2
  DT <: TypeE{6} && return m == 2
  DT <: TypeE{7} && return m == 1
  DT <: TypeE{8} && return m == 8
  DT <: TypeF4 && return m == 4
  DT <: TypeG2 && return m == 1
  return false
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Pretty-printing
# ═══════════════════════════════════════════════════════════════════════════════

function Base.show(io::IO, ::MarkedDynkinType{DT,Marked}) where {DT,Marked}
  print(io, "$(Lie._type_name(DT)) / P_{$(join(Marked, ","))}")
end

function Base.show(io::IO, ::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
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
function marked_dynkin_diagram(::Type{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  # Get the standard diagram and replace ○ at marked positions with ×
  diagram = dynkin_diagram(DT)
  lines = split(diagram, '\n')

  marked_set = Set(Marked)
  R = rank(DT)

  # For simple types, the first line has the nodes
  # We need to replace the i-th ○ with × if i ∈ Marked
  if DT <: SimpleDynkinType
    # Handle D-type specially (has multi-line structure)
    if DT <: TypeD
      return _marked_diagram_D(DT, Marked)
    end
    if DT <: TypeE
      return _marked_diagram_E(DT, Marked)
    end

    # For linear diagrams (A, B, C, F4, G2), first line has the nodes
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

  # For product types, return diagram with marked annotations
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
