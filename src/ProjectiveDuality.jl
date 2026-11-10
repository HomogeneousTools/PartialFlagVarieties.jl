export projective_dual_degree, projective_dual_codimension, projective_dual_defect

struct _ProjectiveDualData
  degree::Union{Nothing,BigInt}
  defect::Union{Nothing,Int}
  section_eulers::Vector{BigInt}
end

const _PROJECTIVE_DUAL_CACHE = LRU{
  Tuple{Union{PartialFlagVariety,ZeroLocus},CompletelyReducibleBundle},
  _ProjectiveDualData,
}(; maxsize=_cache_maxsize(_default_cache_budget(), _DEFAULT_STRUCTURAL_FRAC * 0.05),
  by=Base.summarysize)

# Only compute each section once. Sections beyond the dimension are empty.
function _linear_section_euler!(eulers, X, L, k::Int)
  k > dimension(X) && return BigInt(0)
  E = X isa ZeroLocus ? defining_bundle(X) : zero_bundle(X)
  while length(eulers) <= k
    j = length(eulers)
    push!(eulers, euler_characteristic(zero_locus(E + j * L)))
  end
  eulers[k + 1]
end

function _compute_projective_dual_data(X, L)
  n = dimension(X)
  eulers = BigInt[euler_characteristic(X)]
  for k in 0:n
    e0 = _linear_section_euler!(eulers, X, L, k)
    e1 = _linear_section_euler!(eulers, X, L, k + 1)
    e2 = _linear_section_euler!(eulers, X, L, k + 2)
    polar_degree = (isodd(n - k) ? -1 : 1) * (e0 - 2 * e1 + e2)
    polar_degree < 0 && throw(
      ArgumentError(
        "A negative polar degree was computed. The source and its general sections must be smooth regular zero loci."
      ),
    )
    iszero(polar_degree) && continue
    # A linear projective space has all earlier polar degrees zero and a_n=1.
    # In its own linear span its dual is empty, including a single point.
    k == n && polar_degree == 1 && return _ProjectiveDualData(nothing, nothing, eulers)
    return _ProjectiveDualData(polar_degree, k, eulers)
  end
  # This also covers an empty regular zero locus.
  _ProjectiveDualData(nothing, nothing, eulers)
end

function _projective_dual_data(
  X::Union{PartialFlagVariety,ZeroLocus}, L::CompletelyReducibleBundle
)
  ambient = X isa ZeroLocus ? ambient_variety(X) : X
  _check_polarization(ambient, L)
  is_ample_line_bundle(L) || throw(ArgumentError("The ambient polarization must be ample."))
  # Geometric equality ignores the stored polarization, so include L in the key.
  get!(_PROJECTIVE_DUAL_CACHE, (X, L)) do
    _compute_projective_dual_data(X, L)
  end
end

"""
    projective_dual_degree(X[, L]) -> Union{BigInt,Nothing}

Return the degree of the projective dual of a smooth partial flag variety or
smooth regular zero locus `X`, including when the dual has positive defect.
Omitting `L` uses [`polarization`](@ref) of `X`. An explicit `L` must be an
ample ambient line bundle, or its recorded restriction on a zero locus.
Its restriction specifies a very ample embedding over the complex numbers.

The computation uses **topological** Euler characteristics of general linear
sections, via the classical class formula and the jet-bundle criterion in
Tevelev's *Projective Duality and Homogeneous Spaces*, Theorem 7.1 and
Example 7.5 ([arXiv numbering](https://arxiv.org/abs/math/0112028)).
Smoothness and regularity are assumed, as for [`zero_locus`](@ref).
The three dual invariants share a private cache, including the section Euler
characteristics, and stop at the first nonzero polar degree.

All three functions return `nothing` for an empty dual. A linear projective
space, including a point, is considered in its own linear span and has empty
dual. An empty regular zero locus also returns `nothing`.
See [`projective_dual_codimension`](@ref) and [`projective_dual_defect`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> projective_dual_degree(Gr(2, 5))
5

julia> projective_dual_degree(projective_space(2), line_bundle(projective_space(2), 3))
12

julia> X = projective_space(4); Z = zero_locus(line_bundle(X, 5));

julia> projective_dual_degree(Z)
320

julia> projective_dual_degree(projective_space(3)) === nothing
true
```
"""
projective_dual_degree(
  X::Union{PartialFlagVariety,ZeroLocus}, L::CompletelyReducibleBundle
) =
  _projective_dual_data(X, L).degree

"""
    projective_dual_codimension(X[, L]) -> Union{Int,Nothing}

Return the codimension of the projective dual in the dual projective space,
using [`polarization`](@ref) of `X` when `L` is omitted.
For a nonempty dual this is [`projective_dual_defect`](@ref) plus one.
The assumptions, allowed polarizations, cache, and empty-dual convention are
those of [`projective_dual_degree`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> projective_dual_codimension(Gr(2, 5))
3

julia> projective_dual_codimension(quadric(3))
1
```
"""
function projective_dual_codimension(
  X::Union{PartialFlagVariety,ZeroLocus}, L::CompletelyReducibleBundle
)
  defect = _projective_dual_data(X, L).defect
  defect === nothing ? nothing : defect + 1
end

"""
    projective_dual_defect(X[, L]) -> Union{Int,Nothing}

Return the dual defect: the codimension of the nonempty projective dual
minus one. A hypersurface dual has defect zero.
Omitting `L` uses [`polarization`](@ref) of `X`.
The assumptions, allowed polarizations, cache, and empty-dual convention are
those of [`projective_dual_degree`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> projective_dual_defect(OGr(5, 10))
4
```
"""
projective_dual_defect(
  X::Union{PartialFlagVariety,ZeroLocus}, L::CompletelyReducibleBundle
) =
  _projective_dual_data(X, L).defect

projective_dual_degree(X::Union{PartialFlagVariety,ZeroLocus}) =
  projective_dual_degree(X, _ambient_polarization(X))
projective_dual_codimension(X::Union{PartialFlagVariety,ZeroLocus}) =
  projective_dual_codimension(X, _ambient_polarization(X))
projective_dual_defect(X::Union{PartialFlagVariety,ZeroLocus}) =
  projective_dual_defect(X, _ambient_polarization(X))

projective_dual_degree(Z::ZeroLocus, L::ZeroLocusBundle) =
  projective_dual_degree(Z, _ambient_line_bundle(Z, L))
projective_dual_codimension(Z::ZeroLocus, L::ZeroLocusBundle) =
  projective_dual_codimension(Z, _ambient_line_bundle(Z, L))
projective_dual_defect(Z::ZeroLocus, L::ZeroLocusBundle) =
  projective_dual_defect(Z, _ambient_line_bundle(Z, L))
