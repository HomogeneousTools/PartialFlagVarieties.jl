# ═══════════════════════════════════════════════════════════════════════════════
#  PartialFlagVariety — the main user-facing type for partial flag varieties
#
#  A PartialFlagVariety{MDT} is a thin wrapper around a MarkedDynkinType,
#  providing the primary API for dimension, Picard rank, Betti numbers, etc.
#  Mathematical computations are implemented at the type level via @generated
#  functions, leveraging the compile-time type parameters of MDT.
# ═══════════════════════════════════════════════════════════════════════════════

export PartialFlagVariety
export partial_flag_variety, full_flag_variety
export dynkin_type, dimension, picard_rank
export euler_characteristic, betti_numbers
export is_generalized_grassmannian, is_cominuscule, is_minuscule
export is_adjoint, is_coadjoint, is_full_flag

"""
    PartialFlagVariety{MDT}

A partial flag variety ``G/P`` encoded by a [`MarkedDynkinType`](@ref) `MDT`,
with an optional human-readable name for display.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = PartialFlagVariety{MarkedDynkinType{TypeA{3}, (2,)}}("Gr(2, 4)")
Gr(2, 4)

julia> dimension(X)
4

julia> betti_numbers(X)
5-element Vector{BigInt}:
 1
 1
 2
 1
 1
```
"""
struct PartialFlagVariety{MDT<:MarkedDynkinType}
  name::String
end

PartialFlagVariety{MDT}() where {MDT} = PartialFlagVariety{MDT}("")

# ═══════════════════════════════════════════════════════════════════════════════
#  Constructors
# ═══════════════════════════════════════════════════════════════════════════════

"""
    partial_flag_variety(::Type{DT}, marked_nodes, name="") -> PartialFlagVariety

Construct a partial flag variety from a Dynkin type and a set of marked
(nonparabolic / crossed-out) nodes.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> partial_flag_variety(TypeA{3}, (2,))
A3 / P_{2}

julia> partial_flag_variety(TypeB{4}, (1, 3), "custom name")
custom name
```
"""
function partial_flag_variety(::Type{DT}, marked::NTuple{K,Int}, name::String="") where {DT<:DynkinType,K}
  sorted = Tuple(sort(collect(marked)))
  MDT = MarkedDynkinType{DT,sorted}
  return PartialFlagVariety{MDT}(name)
end

function partial_flag_variety(::Type{DT}, marked::Vector{Int}, name::String="") where {DT<:DynkinType}
  return partial_flag_variety(DT, Tuple(sort(marked)), name)
end

function partial_flag_variety(::Type{DT}, marked::Int, name::String="") where {DT<:DynkinType}
  return partial_flag_variety(DT, (marked,), name)
end

"""
    PartialFlagVariety(dynkin_str::AbstractString, marked::Vector{<:Integer}) -> PartialFlagVariety

Construct a partial flag variety from a Dynkin type string and a set of marked nodes.

The string is parsed by [`parse_dynkin_type`](@ref) (e.g., `"A3"`, `"A2xB3"`, `"E6"`).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = PartialFlagVariety("A3", [2])
A3 / P_{2}

julia> dimension(X)
4

julia> Y = PartialFlagVariety("A2xB3", [1, 3, 4])
A2 × B3 / P_{1,3,4}

julia> picard_rank(Y)
3
```
"""
function PartialFlagVariety(s::AbstractString, marked::Vector{<:Integer})
  DT = parse_dynkin_type(s)
  return partial_flag_variety(DT, Vector{Int}(marked))
end

"""
    full_flag_variety(::Type{DT}, name="") -> PartialFlagVariety

Construct the full (complete) flag variety ``G/B``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = full_flag_variety(TypeA{2});

julia> is_full_flag(X)
true

julia> dimension(X)
3
```
"""
function full_flag_variety(::Type{DT}, name::String="") where {DT<:DynkinType}
  R = rank(DT)
  marked = Tuple(1:R)
  return partial_flag_variety(DT, marked, name)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Display
# ═══════════════════════════════════════════════════════════════════════════════

function Base.show(io::IO, X::PartialFlagVariety{MDT}) where {MDT}
  if X.name != ""
    print(io, X.name)
  else
    print(io, MDT())
  end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type-level accessors (delegated to MDT)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    marked_type(::PartialFlagVariety{MDT}) -> Type{MDT}

Return the `MarkedDynkinType` type parameter.
"""
marked_type(::PartialFlagVariety{MDT}) where {MDT} = MDT

dynkin_type(::PartialFlagVariety{MDT}) where {MDT} = _ambient_type(MDT)
marked_nodes(X::PartialFlagVariety) = marked_nodes(marked_type(X))
unmarked_nodes(X::PartialFlagVariety) = unmarked_nodes(marked_type(X))
Lie.rank(X::PartialFlagVariety) = rank(marked_type(X))

# ═══════════════════════════════════════════════════════════════════════════════
#  Dimension and topological invariants (primary methods)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dimension(::PartialFlagVariety) -> Int

Return the dimension of the partial flag variety ``G/P``.

This equals the number of positive roots of ``G`` that are not roots of the Levi:
``\\dim(G/P) = |\\Phi^+_G| - |\\Phi^+_L|``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> dimension(Gr(2, 4))
4

julia> dimension(partial_flag_variety(TypeA{4}, (1,)))
4

julia> dimension(partial_flag_variety(TypeE{6}, (1,)))
16
```
"""
@generated function dimension(::PartialFlagVariety{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  unmarked = [i for i in 1:R if !(i in Marked)]
  n_pos_G = n_positive_roots(DT)
  if isempty(unmarked)
    n_pos_L = 0
  else
    C = Lie._cartan_matrix_data(DT)
    C_sub = C[unmarked, unmarked]
    ct = cartan_type(C_sub)
    lt = _cartan_type_to_dynkin_type(ct)
    n_pos_L = lt === nothing ? 0 : n_positive_roots(lt)
  end
  d = n_pos_G - n_pos_L
  return :($d)
end

"""
    picard_rank(::PartialFlagVariety) -> Int

Return the Picard rank of ``G/P``, which equals the number of marked nodes.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> picard_rank(Gr(2, 4))
1

julia> picard_rank(partial_flag_variety(TypeA{3}, (1, 3)))
2
```
"""
picard_rank(X::PartialFlagVariety) = central_rank(marked_type(X))

"""
    euler_characteristic(::PartialFlagVariety) -> BigInt

Return the Euler characteristic ``\\chi(G/P) = |W_G| / |W_L|``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> euler_characteristic(Gr(2, 4))
6

julia> euler_characteristic(partial_flag_variety(TypeA{4}, (1,)))
5

julia> euler_characteristic(partial_flag_variety(TypeE{6}, (1,)))
27
```
"""
@generated function euler_characteristic(::PartialFlagVariety{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  unmarked = [i for i in 1:R if !(i in Marked)]
  wG = weyl_order(DT)
  if isempty(unmarked)
    wL = BigInt(1)
  else
    C = Lie._cartan_matrix_data(DT)
    C_sub = C[unmarked, unmarked]
    ct = cartan_type(C_sub)
    lt = _cartan_type_to_dynkin_type(ct)
    wL = lt === nothing ? BigInt(1) : weyl_order(lt)
  end
  chi = wG ÷ wL
  return :($chi)
end

"""
    betti_numbers(::PartialFlagVariety) -> Vector{BigInt}

Compute the Betti numbers of ``G/P`` using the Poincaré polynomial formula:

``P(t) = \\prod_i \\frac{1 - t^{d_i^G}}{1 - t} \\bigg/ \\prod_j \\frac{1 - t^{d_j^L}}{1 - t}``

where ``d_i^G`` and ``d_j^L`` are the degrees of the fundamental invariants of
the Weyl groups of ``G`` and ``L``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> betti_numbers(partial_flag_variety(TypeA{2}, (1,)))
3-element Vector{BigInt}:
 1
 1
 1

julia> betti_numbers(Gr(2, 4))
5-element Vector{BigInt}:
 1
 1
 2
 1
 1
```
"""
@generated function betti_numbers(::PartialFlagVariety{MarkedDynkinType{DT,Marked}}) where {DT,Marked}
  R = rank(DT)
  unmarked = [i for i in 1:R if !(i in Marked)]

  degs_G = collect(degrees_fundamental_invariants(DT))

  if isempty(unmarked)
    degs_L = Int[]
  else
    C = Lie._cartan_matrix_data(DT)
    C_sub = C[unmarked, unmarked]
    ct = cartan_type(C_sub)
    lt = _cartan_type_to_dynkin_type(ct)
    degs_L = lt === nothing ? Int[] : collect(degrees_fundamental_invariants(lt))
  end

  # Compute Poincaré polynomial as a coefficient vector
  function expand_factor(d)
    return ones(BigInt, d)
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
#  Levi and structural queries
# ═══════════════════════════════════════════════════════════════════════════════

levi_type(X::PartialFlagVariety) = levi_type(marked_type(X))
levi_rank(X::PartialFlagVariety) = levi_rank(marked_type(X))
central_rank(X::PartialFlagVariety) = central_rank(marked_type(X))

# ═══════════════════════════════════════════════════════════════════════════════
#  Classification predicates
# ═══════════════════════════════════════════════════════════════════════════════

is_generalized_grassmannian(X::PartialFlagVariety) = is_generalized_grassmannian(marked_type(X))
is_cominuscule(X::PartialFlagVariety) = is_cominuscule(marked_type(X))
is_minuscule(X::PartialFlagVariety) = is_minuscule(marked_type(X))
is_adjoint(X::PartialFlagVariety) = is_adjoint(marked_type(X))
is_coadjoint(X::PartialFlagVariety) = is_coadjoint(marked_type(X))
is_full_flag(X::PartialFlagVariety) = is_full_flag(marked_type(X))

# ═══════════════════════════════════════════════════════════════════════════════
#  Diagram
# ═══════════════════════════════════════════════════════════════════════════════

marked_dynkin_diagram(X::PartialFlagVariety) = marked_dynkin_diagram(marked_type(X))

# ═══════════════════════════════════════════════════════════════════════════════
#  Equality
# ═══════════════════════════════════════════════════════════════════════════════

Base.:(==)(X₁::PartialFlagVariety{MDT}, X₂::PartialFlagVariety{MDT}) where {MDT} = true
Base.:(==)(::PartialFlagVariety, ::PartialFlagVariety) = false
