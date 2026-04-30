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

User-facing wrapper for a partial flag variety ``G/P``, storing its runtime
[`MarkedDynkinType`](@ref) together with an optional display name.
"""
struct PartialFlagVariety
  name::String
  mdt::MarkedDynkinType
end

PartialFlagVariety(mdt::MarkedDynkinType, name::String="") = PartialFlagVariety(name, mdt)

"""
    partial_flag_variety(DT::Type{<:DynkinType}, marked, name="") -> PartialFlagVariety
    partial_flag_variety(s::AbstractString, marked, name="") -> PartialFlagVariety

Construct the partial flag variety ``G/P_I`` with marked nodes `marked`.

Pass either a Dynkin type `DT` such as `TypeA{4}` or a Dynkin-type string such
as `"A4"` or `"A2xB3"`. The marked nodes may be given as a single integer, a
tuple, or a vector; vector input is sorted and converted to `Int`. The optional
`name` is only used for display.

For a convenience alias, `PartialFlagVariety("A3", [2])` is equivalent to
`partial_flag_variety("A3", [2])`.

!!! note
    The one-argument constructor `PartialFlagVariety("31")` is different: it
    decodes a ZeroLocus62 label rather than parsing a Dynkin type string.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (2,));

julia> dimension(X)
6

julia> marked_nodes(partial_flag_variety("A4", 2))
(2,)

julia> marked_nodes(partial_flag_variety("A3", [1, 3]))
(1, 3)

julia> marked_nodes(PartialFlagVariety("A3", [2]))
(2,)
```
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

"""
    PartialFlagVariety(s::AbstractString, marked::Vector{<:Integer}) -> PartialFlagVariety

Convenience alias for [`partial_flag_variety(s, marked)`](@ref).

This two-argument form parses `s` as a Dynkin-type string such as `"A3"` or
`"A2xB3"`. It is distinct from the one-argument constructor
[`PartialFlagVariety(label::AbstractString)`](@ref), which decodes a
ZeroLocus62 label.
"""
PartialFlagVariety(s::AbstractString, marked::Vector{<:Integer}) = partial_flag_variety(
  parse_dynkin_type(s), Vector{Int}(marked)
)

"""Construct the full flag variety ``G/B`` of type `DT`."""
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

"""
    marked_type(X::PartialFlagVariety) -> MarkedDynkinType
    marked_dynkin_type(X::PartialFlagVariety) -> MarkedDynkinType

Return the runtime [`MarkedDynkinType`](@ref) attached to the variety `X`.
`marked_type` is the short alias.
"""
marked_type(X::PartialFlagVariety) = X.mdt
marked_dynkin_type(X::PartialFlagVariety) = X.mdt
dynkin_type(X::PartialFlagVariety) = dynkin_type(X.mdt)
marked_nodes(X::PartialFlagVariety) = marked_nodes(X.mdt)
unmarked_nodes(X::PartialFlagVariety) = unmarked_nodes(X.mdt)
Lie.rank(X::PartialFlagVariety) = rank(dynkin_type(X))

"""
    dimension(X::PartialFlagVariety) -> Int

Return ``\\dim(G/P) = |\\Phi_G^+| - |\\Phi_L^+|``.
"""
dimension(X::PartialFlagVariety) = dimension(marked_dynkin_type(X))

"""
    picard_rank(X::PartialFlagVariety) -> Int

Return the Picard rank of `X`, equal to the number of marked nodes
``|I|``.
"""
picard_rank(X::PartialFlagVariety) = central_rank(marked_dynkin_type(X))

levi_type(X::PartialFlagVariety) = levi_type(marked_dynkin_type(X))
levi_rank(X::PartialFlagVariety) = levi_rank(marked_dynkin_type(X))
central_rank(X::PartialFlagVariety) = central_rank(marked_dynkin_type(X))

"""
    euler_characteristic(X::PartialFlagVariety) -> BigInt

Return the topological Euler characteristic ``\\chi(G/P) = |W_G|/|W_L|``.
"""
function euler_characteristic(X::PartialFlagVariety)
  wG = weyl_order(dynkin_type(X))
  wL = is_borel(marked_dynkin_type(X)) ? BigInt(1) : weyl_order(levi_type(X))
  wG ÷ wL
end

"""
    betti_numbers(X::PartialFlagVariety) -> Vector{BigInt}

Return the Betti numbers of `X` in even degrees, i.e.
``[b_0, b_2, b_4, \\ldots, b_{2d}]``.

Since ``G/P`` has no odd cohomology this completely describes ``\\mathrm{H}^*(X, \\mathbb{Z})``.
Computed from the ratio of the Poincaré polynomials of the Weyl groups
``W_G`` and ``W_L``.
"""
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

"""
    is_generalized_grassmannian(X::PartialFlagVariety) -> Bool

Return `true` for generalized Grassmannians, i.e. varieties ``G/P_i``
with exactly one marked node.
"""
is_generalized_grassmannian(X::PartialFlagVariety) = length(marked_nodes(X)) == 1

"""
    is_full_flag_variety(X::PartialFlagVariety) -> Bool

Return `true` for full flag varieties ``G/B`` (all nodes marked).
"""
is_full_flag_variety(X::PartialFlagVariety) = is_borel(marked_dynkin_type(X))

"""
    is_cominuscule(X::PartialFlagVariety) -> Bool

Return `true` if `X` is a cominuscule flag variety, i.e., a generalized
Grassmannian ``G/P_i`` whose highest root has coefficient 1 at node ``i``.

Cominuscule varieties include all Grassmannians (type A), spinor varieties
``\\mathrm{OGr}(n, 2n)``, Lagrangian Grassmannians ``\\mathrm{LGr}(n, 2n)``,
the two connected components of the orthogonal Grassmannian in type D,
the Cayley plane (``E_6/P_1``, ``E_6/P_6``), and the Freudenthal variety
(``E_7/P_7``).
"""
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

"""
    is_minuscule(X::PartialFlagVariety) -> Bool

Return `true` if `X` is a minuscule flag variety, i.e., a generalized
Grassmannian ``G/P_i`` whose fundamental weight ``\\omega_i`` is minuscule
(all coroot pairings are 0 or 1).

Minuscule varieties are the same as cominuscule for simply laced types, but
differ for types B and C.
"""
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

"""
    is_adjoint(X::PartialFlagVariety) -> Bool

Return `true` if `X` is the adjoint variety of its Lie type, i.e., the
orbit of the highest root line in the projectivisation of the adjoint
representation.
"""
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

"""
    is_coadjoint(X::PartialFlagVariety) -> Bool

Return `true` if `X` is the coadjoint variety of its Lie type, i.e., the
orbit of the highest weight line in the dual of the adjoint representation.
For simply laced types the adjoint and coadjoint coincide.
"""
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

"""
    anticanonical_degrees(X::PartialFlagVariety) -> Vector{Int}

Return the coefficients ``(a_1, \\ldots, a_r)`` of the anticanonical
divisor ``-\\mathrm{K}_X = \\sum_j a_j [D_j]`` in the marked-node basis of
``\\operatorname{Pic}(X)``.

Computed via the formula
``a_j = \\langle 2(\\rho_G - \\rho_P),\\, \\alpha_j^\\vee \\rangle``
where ``\\rho_G`` is the Weyl vector of ``G`` and ``\\rho_P`` is the Weyl
vector of the Levi factor.
"""
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
