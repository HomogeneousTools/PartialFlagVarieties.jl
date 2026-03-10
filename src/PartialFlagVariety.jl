export PartialFlagVariety
export partial_flag_variety, full_flag_variety
export dynkin_type, dimension, picard_rank
export euler_characteristic, betti_numbers
export is_generalized_grassmannian, is_cominuscule, is_minuscule
export is_adjoint, is_coadjoint, is_full_flag_variety
export anticanonical_degrees
export marked_type, marked_dynkin_type, marked_nodes

"""
    PartialFlagVariety(mdt::MarkedDynkinType, name="")

User-facing wrapper for a partial flag variety `G/P`, storing its runtime
[`MarkedDynkinType`](@ref) together with an optional display name.
"""
struct PartialFlagVariety
  name::String
  mdt::MarkedDynkinType
end

PartialFlagVariety(mdt::MarkedDynkinType, name::String="") = PartialFlagVariety(name, mdt)

"""
    partial_flag_variety(DT::Type{<:DynkinType}, marked, name="") -> PartialFlagVariety

Construct the partial flag variety of type `DT/P_marked`.
"""
function partial_flag_variety(
  ::Type{DT}, marked::Tuple, name::String=""
) where {DT<:DynkinType}
  PartialFlagVariety(MarkedDynkinType(DT, marked), name)
end

partial_flag_variety(::Type{DT}, marked::Vector{<:Integer}, name::String="") where {DT<:DynkinType} = partial_flag_variety(
  DT, Tuple(sort(Int.(marked))), name
)

partial_flag_variety(::Type{DT}, marked::Integer, name::String="") where {DT<:DynkinType} = partial_flag_variety(
  DT, (Int(marked),), name
)

function partial_flag_variety(
  s::AbstractString, marked::AbstractVector{<:Integer}, name::String=""
)
  partial_flag_variety(parse_dynkin_type(s), Vector{Int}(marked), name)
end

function partial_flag_variety(s::AbstractString, marked::Integer, name::String="")
  partial_flag_variety(parse_dynkin_type(s), Int(marked), name)
end

PartialFlagVariety(s::AbstractString, marked::Vector{<:Integer}) = partial_flag_variety(
  parse_dynkin_type(s), Vector{Int}(marked)
)

"""Construct the full flag variety `G/B` of type `DT`."""
function full_flag_variety(::Type{DT}, name::String="") where {DT<:DynkinType}
  partial_flag_variety(DT, Tuple(1:rank(DT)), name)
end

function Base.show(io::IO, X::PartialFlagVariety)
  if X.name != ""
    print(io, X.name)
  else
    print(io, marked_dynkin_type(X))
  end
end

"""Return the runtime marked Dynkin type attached to `X`."""
marked_type(X::PartialFlagVariety) = X.mdt
marked_dynkin_type(X::PartialFlagVariety) = X.mdt
dynkin_type(X::PartialFlagVariety) = dynkin_type(X.mdt)
marked_nodes(X::PartialFlagVariety) = marked_nodes(X.mdt)
unmarked_nodes(X::PartialFlagVariety) = unmarked_nodes(X.mdt)
Lie.rank(X::PartialFlagVariety) = rank(dynkin_type(X))

"""Return the dimension of `X`."""
dimension(X::PartialFlagVariety) = dimension(marked_dynkin_type(X))

"""Return the Picard rank of `X`."""
picard_rank(X::PartialFlagVariety) = central_rank(marked_dynkin_type(X))

levi_type(X::PartialFlagVariety) = levi_type(marked_dynkin_type(X))
levi_rank(X::PartialFlagVariety) = levi_rank(marked_dynkin_type(X))
central_rank(X::PartialFlagVariety) = central_rank(marked_dynkin_type(X))

"""Return the topological Euler characteristic of `X`."""
function euler_characteristic(X::PartialFlagVariety)
  wG = weyl_order(dynkin_type(X))
  wL = is_borel(marked_dynkin_type(X)) ? BigInt(1) : weyl_order(levi_type(X))
  wG ÷ wL
end

"""Return the Betti numbers of `X` in even degrees."""
function betti_numbers(X::PartialFlagVariety)
  degs_G = collect(degrees_fundamental_invariants(dynkin_type(X)))
  degs_L = if is_borel(marked_dynkin_type(X))
    Int[]
  else
    collect(degrees_fundamental_invariants(levi_type(X)))
  end

  num = BigInt[1]
  for d in degs_G
    num = _poly_mul(num, ones(BigInt, d))
  end

  den = BigInt[1]
  for d in degs_L
    den = _poly_mul(den, ones(BigInt, d))
  end

  _poly_div(num, den)
end

function _poly_mul(a, b)
  la, lb = length(a), length(b)
  result = zeros(BigInt, la + lb - 1)
  for i in 1:la, j in 1:lb
    result[i + j - 1] += a[i] * b[j]
  end
  result
end

function _poly_div(a, b)
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
  result
end

"""Return `true` for generalized Grassmannians, i.e. varieties with one marked node."""
is_generalized_grassmannian(X::PartialFlagVariety) = length(marked_nodes(X)) == 1

"""Return `true` for full flag varieties `G/B`."""
is_full_flag_variety(X::PartialFlagVariety) = is_borel(marked_dynkin_type(X))

function is_cominuscule(X::PartialFlagVariety)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  DT <: SimpleDynkinType || return false
  length(marked) != 1 && return false
  m = marked[1]
  R = rank(X)
  DT <: TypeA && return true
  DT <: TypeB && return m == 1
  DT <: TypeC && return m == R
  DT <: TypeD && return (m == 1 || m == R - 1 || m == R)
  DT <: TypeE{6} && return (m == 1 || m == 6)
  DT <: TypeE{7} && return m == 7
  false
end

function is_minuscule(X::PartialFlagVariety)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  DT <: SimpleDynkinType || return false
  length(marked) != 1 && return false
  m = marked[1]
  R = rank(X)
  DT <: TypeA && return true
  DT <: TypeB && return m == R
  DT <: TypeC && return m == 1
  DT <: TypeD && return (m == 1 || m == R - 1 || m == R)
  DT <: TypeE{6} && return (m == 1 || m == 6)
  DT <: TypeE{7} && return m == 7
  false
end

function is_adjoint(X::PartialFlagVariety)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  R = rank(X)
  DT <: SimpleDynkinType || return false
  if DT <: TypeA
    return length(marked) == 2 && marked == (1, R)
  end
  length(marked) != 1 && return false
  m = marked[1]
  DT <: TypeB && return m == 2
  DT <: TypeC && return m == 1
  DT <: TypeD && R >= 4 && return m == 2
  DT <: TypeE{6} && return m == 2
  DT <: TypeE{7} && return m == 1
  DT <: TypeE{8} && return m == 8
  DT <: TypeF4 && return m == 1
  DT <: TypeG2 && return m == 2
  false
end

function is_coadjoint(X::PartialFlagVariety)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  R = rank(X)
  DT <: SimpleDynkinType || return false
  if DT <: TypeA
    return length(marked) == 2 && marked == (1, R)
  end
  length(marked) != 1 && return false
  m = marked[1]
  DT <: TypeB && return m == 1
  DT <: TypeC && return m == 2
  DT <: TypeD && R >= 4 && return m == 2
  DT <: TypeE{6} && return m == 2
  DT <: TypeE{7} && return m == 1
  DT <: TypeE{8} && return m == 8
  DT <: TypeF4 && return m == 4
  DT <: TypeG2 && return m == 1
  false
end

marked_dynkin_diagram(X::PartialFlagVariety) = marked_dynkin_diagram(marked_dynkin_type(X))

"""Return the coefficients of `-K_X` in the marked-node basis of `Pic(X)`."""
function anticanonical_degrees(X::PartialFlagVariety)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  unmarked = unmarked_nodes(X)
  n = length(unmarked)
  C = Lie._cartan_matrix_data(DT)

  if n == 0
    return fill(2, length(marked))
  end

  C_L = Rational{Int}[C[unmarked[p], unmarked[q]] for p in 1:n, q in 1:n]
  x_L = C_L \ ones(Rational{Int}, n)

  Int[round(Int, 2 - 2 * sum(C[i, unmarked[q]] * x_L[q] for q in 1:n)) for i in marked]
end

Base.:(==)(X₁::PartialFlagVariety, X₂::PartialFlagVariety) =
  marked_dynkin_type(X₁) == marked_dynkin_type(X₂)
Base.hash(X::PartialFlagVariety, h::UInt) = hash(marked_dynkin_type(X), h)
