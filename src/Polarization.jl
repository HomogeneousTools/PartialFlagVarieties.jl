export polarization, set_polarization!

# Store Picard degrees rather than a bundle referring back to its own base.
# A zero locus without an override follows its ambient variety's polarization.
_polarization_degrees(X::PartialFlagVariety) =
  something(X.polarization, ntuple(_ -> 1, picard_rank(X)))
_polarization_degrees(Z::ZeroLocus) =
  something(Z.polarization, _polarization_degrees(ambient_variety(Z)))

_ambient_polarization(X::PartialFlagVariety) = polarization(X)
_ambient_polarization(Z::ZeroLocus) =
  line_bundle(ambient_variety(Z), collect(_polarization_degrees(Z)))

"""
    polarization(X::Union{PartialFlagVariety,ZeroLocus})

Return the chosen ample line bundle on `X`. Unless set explicitly with
[`set_polarization!`](@ref), a partial flag variety uses degree one at every
marked node, and a zero locus inherits the current ambient polarization.
The result lives on `X`, including when `X` is a zero locus.

Hilbert polynomials use this choice when their line-bundle argument is omitted.
Products carry the external tensor product of the factor polarizations.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(2); set_polarization!(X, 2);

julia> picard_degrees(polarization(X))
1-element Vector{Int64}:
 2

julia> hilbert_polynomial(structure_sheaf(X))
Polynomial(1//1 + 3//1*t + 2//1*t^2)
```
"""
polarization(X::PartialFlagVariety) = line_bundle(X, collect(_polarization_degrees(X)))
polarization(Z::ZeroLocus) = restrict(Z, _ambient_polarization(Z))

# Recognize an explicitly recorded restriction of an ambient line bundle.
# Arbitrary intrinsic presentations need not extend to the ambient variety.
function _ambient_line_bundle(Z::ZeroLocus, L::CompletelyReducibleBundle)
  _check_polarization(ambient_variety(Z), L)
  L
end

function _ambient_line_bundle(Z::ZeroLocus, L::ZeroLocusBundle)
  variety(L) == Z || throw(ArgumentError("The line bundle must live on this zero locus."))
  terms = L.presentation.terms
  if length(terms) == 1 && haskey(terms, 0) && length(terms[0]) == 1
    ambient = only(terms[0])
    ambient isa CompletelyReducibleBundle && return _ambient_line_bundle(Z, ambient)
  end
  throw(
    ArgumentError(
      "Supply a line bundle recorded as a restriction from the ambient variety."
    ),
  )
end

"""
    set_polarization!(X, L)
    set_polarization!(X, degrees)
    set_polarization!(X, nothing)

Choose an ample line bundle for a [`PartialFlagVariety`](@ref) or
[`ZeroLocus`](@ref), returning `X`. `degrees` is a vector of positive Picard
degrees, or an integer for an ambient variety of Picard rank one.
For a zero locus, `L` may be an ambient line bundle or its recorded restriction.
Ampleness is checked on the ambient variety.

Passing `nothing` removes the override: a partial flag variety returns to
degree one at every marked node, and a zero locus follows its ambient variety.
The constructors `partial_flag_variety(...; polarization=degrees)` and
`zero_locus(E; polarization=L)` also accept these choices.

Equality, hashing, and geometric labels ignore this mutable choice; two
presentations of the same variety can have different polarizations.
An explicit line-bundle argument to `hilbert_polynomial` overrides the stored
choice for that call without changing it.
"""
function set_polarization!(X::PartialFlagVariety, L::CompletelyReducibleBundle)
  _check_polarization(X, L)
  is_ample_line_bundle(L) || throw(ArgumentError("The polarization must be ample."))
  X.polarization = Tuple(picard_degrees(L))
  X
end

function set_polarization!(
  Z::ZeroLocus, L::Union{CompletelyReducibleBundle,ZeroLocusBundle}
)
  ambient = _ambient_line_bundle(Z, L)
  is_ample_line_bundle(ambient) ||
    throw(ArgumentError("The ambient polarization must be ample."))
  Z.polarization = Tuple(picard_degrees(ambient))
  Z
end

function set_polarization!(
  X::Union{PartialFlagVariety,ZeroLocus}, degrees::Union{Integer,AbstractVector{<:Integer}}
)
  ambient = X isa ZeroLocus ? ambient_variety(X) : X
  d = degrees isa Integer ? degrees : collect(degrees)
  set_polarization!(X, line_bundle(ambient, d))
end

function set_polarization!(X::Union{PartialFlagVariety,ZeroLocus}, ::Nothing)
  X.polarization = nothing
  X
end

hilbert_polynomial(Z::ZeroLocus, L::ZeroLocusBundle) =
  hilbert_polynomial(Z, _ambient_line_bundle(Z, L))
hilbert_polynomial(F::ZeroLocusBundle, L::ZeroLocusBundle) =
  hilbert_polynomial(F, _ambient_line_bundle(variety(F), L))
