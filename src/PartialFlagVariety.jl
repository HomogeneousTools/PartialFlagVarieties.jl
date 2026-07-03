export PartialFlagVariety
export partial_flag_variety, full_flag_variety, product
export dynkin_type, dimension, picard_rank
export euler_characteristic, betti_numbers
export is_generalized_grassmannian, is_cominuscule, is_minuscule, is_exceptional_type
export is_adjoint, is_coadjoint, is_full_flag_variety
export is_orthogonal_grassmannian, is_projective_space, is_quadric
export anticanonical_degrees
export marked_dynkin_type, marked_nodes

"""
    PartialFlagVariety(mdt::MarkedDynkinType, name="")

User-facing wrapper for a partial flag variety ``G/P``, storing its runtime
[`MarkedDynkinType`](@ref) together with an optional display name.

Most users should create these through the named constructors (`Gr`, `quadric`,
`cayley_plane`, ...) or via [`partial_flag_variety`](@ref), rather than by
assembling a `MarkedDynkinType` manually.
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

The convention is that **marked nodes are the crossed-out / nonparabolic
nodes** of the Dynkin diagram.

For a convenience alias, `PartialFlagVariety("A3", [2])` is equivalent to
`partial_flag_variety("A3", [2])`.

!!! note
    The one-argument constructor `PartialFlagVariety("31")` is different: it
    decodes a ZeroLocus62 label rather than parsing a Dynkin type string.

# Examples
```jldoctest
julia> using PartialFlagVarieties

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

partial_flag_variety(
  ::Type{DT}, marked::Vector{<:Integer}, name::String=""
) where {DT<:DynkinType} = partial_flag_variety(DT, Tuple(marked), name)

partial_flag_variety(::Type{DT}, marked::Integer, name::String="") where {DT<:DynkinType} =
  partial_flag_variety(
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

"""
    full_flag_variety(::Type{DT}, name="") -> PartialFlagVariety

Construct the full flag variety ``G/B`` of type `DT`, i.e. the case where every
simple root is marked.
"""
function full_flag_variety(::Type{DT}, name::String="") where {DT<:DynkinType}
  partial_flag_variety(DT, Tuple(1:rank(DT)), name)
end

"""
    _flatten_dynkin_factors(DT::Type{<:DynkinType}) -> Vector{DataType}

Flatten a simple or nested product Dynkin type into its ordered simple factors.
"""
function _flatten_dynkin_factors(::Type{DT}) where {DT<:SimpleDynkinType}
  DataType[DT]
end

function _flatten_dynkin_factors(::Type{DT}) where {DT<:ProductDynkinType}
  result = DataType[]
  for factor in DT.parameters[1].parameters
    append!(result, _flatten_dynkin_factors(factor))
  end
  result
end

"""
    _combine_dynkin_factors(factors) -> Type{<:DynkinType}

Rebuild a Dynkin type from an ordered list of simple factors.
"""
function _combine_dynkin_factors(factors::AbstractVector{<:DataType})
  isempty(factors) && throw(ArgumentError("Need at least one Dynkin factor."))
  length(factors) == 1 && return first(factors)
  ProductDynkinType{Tuple{factors...}}
end

"""
    _product_marked_dynkin_type(varieties...) -> MarkedDynkinType

Construct the marked Dynkin type of a product ambient by concatenating the
simple factors of each variety and shifting marked nodes by cumulative rank.
"""
function _product_marked_dynkin_type(varieties::Vararg{PartialFlagVariety,N}) where {N}
  factor_types = DataType[]
  marked = Int[]
  offset = 0
  for X in varieties
    append!(factor_types, _flatten_dynkin_factors(dynkin_type(X)))
    append!(marked, (offset + m for m in marked_nodes(X)))
    offset += rank(X)
  end
  DT = _combine_dynkin_factors(factor_types)
  MarkedDynkinType(DT, Tuple(marked))
end

"""
    product(X::PartialFlagVariety, Y::PartialFlagVariety, Zs::PartialFlagVariety...) -> PartialFlagVariety

Construct the product of partial flag varieties.

The ambient Dynkin factors are concatenated in order, and the marked nodes of
later factors are shifted by the cumulative ranks of the earlier factors.
This is also available through the `*` operator.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = product(projective_space(1), projective_space(2));

julia> dimension(X)
3

julia> picard_rank(X)
2
```
"""
function product(
  X::PartialFlagVariety, Y::PartialFlagVariety, Zs::PartialFlagVariety...
)
  PartialFlagVariety(_product_marked_dynkin_type(X, Y, Zs...))
end

Base.:*(X::PartialFlagVariety, Y::PartialFlagVariety) = product(X, Y)

function Base.show(io::IO, X::PartialFlagVariety)
  if X.name != ""
    print(io, X.name)
  else
    print(io, marked_dynkin_type(X))
  end
end

"""
    marked_dynkin_type(X::PartialFlagVariety) -> MarkedDynkinType

Return the runtime [`MarkedDynkinType`](@ref) attached to the variety `X`.
"""
marked_dynkin_type(X::PartialFlagVariety) = X.mdt
dynkin_type(X::PartialFlagVariety) = dynkin_type(X.mdt)
marked_nodes(X::PartialFlagVariety) = marked_nodes(X.mdt)
unmarked_nodes(X::PartialFlagVariety) = unmarked_nodes(X.mdt)
rank(X::PartialFlagVariety) = rank(dynkin_type(X))

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

The entry at index `p + 1` is ``b_{2p}``.
"""
function betti_numbers(X::PartialFlagVariety)
  degs_G = collect(degrees_fundamental_invariants(dynkin_type(X)))
  degs_L = if is_borel(marked_dynkin_type(X))
    Int[]
  else
    collect(degrees_fundamental_invariants(levi_type(X)))
  end

  num = foldl((p, d) -> _poly_mul(p, ones(BigInt, d)), degs_G; init=BigInt[1])
  den = foldl((p, d) -> _poly_mul(p, ones(BigInt, d)), degs_L; init=BigInt[1])
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
    is_orthogonal_grassmannian(X::PartialFlagVariety) -> Bool

Return `true` if `X` is a one-marked orthogonal Grassmannian, i.e.
``\\mathrm{B}_n/P_k = \\mathrm{OGr}(k, 2n+1)`` or
``\\mathrm{D}_n/P_k = \\mathrm{OGr}(k, 2n)``.

The two-marked ``\\mathrm{D}_n/P_{n-1, n}`` (the ``(n-1)``-isotropic
Grassmannian, of Picard rank 2) is **not** an orthogonal Grassmannian by this
predicate; [`spinor_bundle`](@ref) accepts it as a separate special case.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_orthogonal_grassmannian(OGr(3, 10))
true

julia> is_orthogonal_grassmannian(partial_flag_variety(TypeD{4}, (3, 4)))
false

julia> is_orthogonal_grassmannian(Gr(2, 5))
false
```
"""
function is_orthogonal_grassmannian(X::PartialFlagVariety)
  DT = dynkin_type(X)
  (DT <: TypeB || DT <: TypeD) || return false
  length(marked_nodes(X)) == 1
end

"""
    is_projective_space(X::PartialFlagVariety) -> Bool

Return `true` if `X` is (isomorphic to) a projective space ``\\mathbb{P}^n``.

In type A these are ``\\mathrm{A}_n/P_1 = \\mathrm{Gr}(1, n+1)`` and its dual
``\\mathrm{A}_n/P_n = \\mathrm{Gr}(n, n+1)``. The remaining cases are low-rank
accidental isomorphisms:

- ``\\mathrm{C}_n/P_1 = \\mathbb{P}^{2n-1}`` — every line is isotropic for a
  symplectic form, so ``\\mathrm{SGr}(1, 2n) = \\mathbb{P}(\\mathbb{C}^{2n})``;
- ``\\mathrm{B}_2/P_2 = \\mathbb{P}^3`` — the spinor variety ``\\mathrm{OGr}(2, 5)``,
  via ``\\mathrm{Spin}(5) = \\mathrm{Sp}(4)``;
- ``\\mathrm{D}_3/P_2 = \\mathrm{D}_3/P_3 = \\mathbb{P}^3`` — the two spinor
  varieties, via ``\\mathrm{D}_3 = \\mathrm{A}_3``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_projective_space(projective_space(4))
true

julia> is_projective_space(SGr(1, 6))   # C₃/P₁ = ℙ⁵
true

julia> is_projective_space(OGr(2, 5))   # B₂/P₂ = ℙ³
true

julia> is_projective_space(Gr(2, 5))
false
```
"""
function is_projective_space(X::PartialFlagVariety)
  DT = dynkin_type(X)
  marked = marked_nodes(X)
  length(marked) != 1 && return false
  m = marked[1]
  R = rank(X)
  DT <: TypeA && return (m == 1 || m == R)
  DT <: TypeC && return m == 1
  DT <: TypeB && return (R == 2 && m == 2)
  DT <: TypeD && return (R == 3 && (m == 2 || m == 3))
  false
end

"""
    is_quadric(X::PartialFlagVariety) -> Bool

Return `true` if `X` is a smooth quadric hypersurface ``Q^n``, i.e.
``\\mathrm{B}_m/P_1`` (odd ``n = 2m-1``) or ``\\mathrm{D}_m/P_1`` (even ``n = 2m-2``).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> is_quadric(quadric(4))
true

julia> is_quadric(OGr(2, 10))
false
```
"""
function is_quadric(X::PartialFlagVariety)
  DT = dynkin_type(X)
  (DT <: TypeB || DT <: TypeD) || return false
  marked_nodes(X) == (1,)
end

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
    is_exceptional_type(X::PartialFlagVariety) -> Bool

Return `true` if the ambient Dynkin type contains an exceptional simple factor
(``\\mathrm{E}_6``, ``\\mathrm{E}_7``, ``\\mathrm{E}_8``, ``\\mathrm{F}_4``, or ``\\mathrm{G}_2``).

For product types, returns `true` if **any** component is exceptional. Note that
``\\mathrm{G}_2/P_1 \\cong Q^5`` geometrically, but is still considered exceptional here since
the Lie type is ``\\mathrm{G}_2``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety("E6", 2);

julia> is_exceptional_type(X)
true
```
"""
is_exceptional_type(X::PartialFlagVariety) = is_exceptional_type(dynkin_type(X))

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

The basis is ordered by [`marked_nodes(X)`](@ref).

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
  C = Semisimple._cartan_matrix_data(DT)

  if n == 0
    return fill(2, length(marked))
  end

  C_L = Rational{Int}[C[unmarked[p], unmarked[q]] for p in 1:n, q in 1:n]
  x_L = C_L \ ones(Rational{Int}, n)

  result = Vector{Int}(undef, length(marked))
  for (p, i) in enumerate(marked)
    total = zero(Rational{Int})
    for q in 1:n
      total += C[i, unmarked[q]] * x_L[q]
    end
    result[p] = Int(2 - 2 * total)
  end
  result
end

Base.:(==)(X₁::PartialFlagVariety, X₂::PartialFlagVariety) =
  marked_dynkin_type(X₁) == marked_dynkin_type(X₂)
Base.hash(X::PartialFlagVariety, h::UInt) = hash(marked_dynkin_type(X), h)
