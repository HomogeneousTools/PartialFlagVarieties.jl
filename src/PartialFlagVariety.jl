# ═══════════════════════════════════════════════════════════════════════════════
#  PartialFlagVariety — named wrapper around MarkedDynkinType
#
#  A PartialFlagVariety{MDT} is a thin wrapper providing a display name
#  and serving as the main user-facing type. All mathematical computations
#  delegate to the type parameter MDT.
# ═══════════════════════════════════════════════════════════════════════════════

export PartialFlagVariety
export partial_flag_variety, full_flag_variety
export dynkin_type, picard_rank
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

julia> V = PartialFlagVariety{MarkedDynkinType{TypeA{3}, (2,)}}("Gr(2, 4)")
Gr(2, 4)

julia> dimension(V)
4

julia> betti_numbers(V)
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
    full_flag_variety(::Type{DT}, name="") -> PartialFlagVariety

Construct the full (complete) flag variety ``G/B``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = full_flag_variety(TypeA{2});

julia> is_full_flag(V)
true

julia> dimension(V)
3
```
"""
function full_flag_variety(::Type{DT}, name::String="") where {DT<:DynkinType}
  R = rank(DT)
  marked = Tuple(1:R)
  return partial_flag_variety(DT, marked, name)
end

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, V::PartialFlagVariety{MDT}) where {MDT}
  if V.name != ""
    print(io, V.name)
  else
    print(io, MDT())
  end
end

# ─── Delegate all methods to MDT ────────────────────────────────────────────

"""
    marked_type(::PartialFlagVariety{MDT}) -> Type{MDT}

Return the `MarkedDynkinType` type parameter.
"""
marked_type(::PartialFlagVariety{MDT}) where {MDT} = MDT

dynkin_type(::PartialFlagVariety{MDT}) where {MDT} = _ambient_type(MDT)
marked_nodes(V::PartialFlagVariety) = marked_nodes(marked_type(V))
unmarked_nodes(V::PartialFlagVariety) = unmarked_nodes(marked_type(V))
Lie.rank(V::PartialFlagVariety) = rank(marked_type(V))

dimension(V::PartialFlagVariety) = dimension(marked_type(V))
picard_rank(V::PartialFlagVariety) = picard_rank(marked_type(V))
euler_characteristic(V::PartialFlagVariety) = euler_characteristic(marked_type(V))
betti_numbers(V::PartialFlagVariety) = betti_numbers(marked_type(V))
levi_type(V::PartialFlagVariety) = levi_type(marked_type(V))
levi_rank(V::PartialFlagVariety) = levi_rank(marked_type(V))
central_rank(V::PartialFlagVariety) = central_rank(marked_type(V))

is_generalized_grassmannian(V::PartialFlagVariety) = is_generalized_grassmannian(marked_type(V))
is_cominuscule(V::PartialFlagVariety) = is_cominuscule(marked_type(V))
is_minuscule(V::PartialFlagVariety) = is_minuscule(marked_type(V))
is_adjoint(V::PartialFlagVariety) = is_adjoint(marked_type(V))
is_coadjoint(V::PartialFlagVariety) = is_coadjoint(marked_type(V))
is_full_flag(V::PartialFlagVariety) = is_full_flag(marked_type(V))

marked_dynkin_diagram(V::PartialFlagVariety) = marked_dynkin_diagram(marked_type(V))

# ─── Equality ────────────────────────────────────────────────────────────────

Base.:(==)(V₁::PartialFlagVariety{MDT}, V₂::PartialFlagVariety{MDT}) where {MDT} = true
Base.:(==)(::PartialFlagVariety, ::PartialFlagVariety) = false

# ─── Product ─────────────────────────────────────────────────────────────────

# TODO: implement product of partial flag varieties via ProductDynkinType
