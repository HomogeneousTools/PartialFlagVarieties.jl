export DegeneracyLocus, DegeneracyLocusBundle, incidence_resolution

"""
    DegeneracyLocus

A general determinantal, Pfaffian, symmetric or orbital degeneracy locus.
Retains the construction data separately from its incidence zero locus.
Use [`incidence_resolution`](@ref) to inspect that presentation.

The constructors assume regularity and expected dimension. When a deeper
stratum is present, cohomology, Hodge and tangent calculations refer to the
incidence resolution, which can differ from its image.
Equality compares construction data and the chosen incidence presentation;
it does not compare sections or recognize isomorphic presentations.
"""
struct DegeneracyLocus{P<:Union{ZeroLocus,RelativeZeroLocus},D<:NamedTuple}
  incidence::P
  construction::Symbol
  data::D
end

"""Return the `ZeroLocus` or `RelativeZeroLocus` presenting a degeneracy locus."""
incidence_resolution(Z::DegeneracyLocus) = Z.incidence
base_variety(Z::DegeneracyLocus) = variety(Z.data.E)
ambient_variety(Z::DegeneracyLocus) = base_variety(Z)
dimension(Z::DegeneracyLocus) = dimension(incidence_resolution(Z))
codimension(Z::DegeneracyLocus) = dimension(ambient_variety(Z)) - dimension(Z)

Base.:(==)(Z::DegeneracyLocus, W::DegeneracyLocus) =
  Z.construction == W.construction && Z.data == W.data && Z.incidence == W.incidence
Base.isequal(Z::DegeneracyLocus, W::DegeneracyLocus) = Z == W
Base.hash(Z::DegeneracyLocus, h::UInt) =
  hash(Z.incidence, hash(Z.data, hash(Z.construction, h)))
Base.show(io::IO, Z::DegeneracyLocus) =
  print(
    io,
    "Degeneracy locus(",
    Z.construction,
    ", dim ",
    dimension(Z),
    " in ",
    ambient_variety(Z),
    ")",
  )

"""A bundle on a degeneracy locus, represented on its incidence resolution."""
struct DegeneracyLocusBundle{
  Z<:DegeneracyLocus,B<:Union{ZeroLocusBundle,RelativeZeroLocusBundle}
} <: Bundle
  variety::Z
  presentation::B

  function DegeneracyLocusBundle(
    Z::DegeneracyLocus, F::Union{ZeroLocusBundle,RelativeZeroLocusBundle}
  )
    variety(F) == incidence_resolution(Z) ||
      throw(ArgumentError("The bundle must live on the incidence resolution."))
    new{typeof(Z),typeof(F)}(Z, F)
  end
end
variety(F::DegeneracyLocusBundle) = F.variety
ambient_variety(F::DegeneracyLocusBundle) = ambient_variety(variety(F))
rank(F::DegeneracyLocusBundle) = rank(F.presentation)
Base.iszero(F::DegeneracyLocusBundle) = iszero(F.presentation)
Base.:(==)(F::DegeneracyLocusBundle, H::DegeneracyLocusBundle) =
  variety(F) == variety(H) && F.presentation == H.presentation
Base.isequal(F::DegeneracyLocusBundle, H::DegeneracyLocusBundle) = F == H
Base.hash(F::DegeneracyLocusBundle, h::UInt) = hash(F.presentation, hash(variety(F), h))
Base.show(io::IO, F::DegeneracyLocusBundle) =
  print(io, "Bundle(rank ", rank(F), " on ", variety(F), ")")

function restrict(Z::DegeneracyLocus, F::Union{_AmbientBundle,RelativeGrassmannianBundle})
  R = incidence_resolution(Z)
  if R isa ZeroLocus && F isa _AmbientBundle && variety(F) == base_variety(Z) &&
    variety(F) != ambient_variety(R)
    F = _lift_bundle_to_product(ambient_variety(R), F, 0)
  end
  DegeneracyLocusBundle(Z, restrict(R, F))
end

line_bundle(Z::DegeneracyLocus, d::Integer) = restrict(Z, line_bundle(base_variety(Z), d))
line_bundle(Z::DegeneracyLocus, degrees::Vector{<:Integer}) =
  restrict(Z, line_bundle(base_variety(Z), degrees))
O(Z::DegeneracyLocus, d::Integer=0) = iszero(d) ? structure_sheaf(Z) : line_bundle(Z, d)
O(Z::DegeneracyLocus, degrees::Vector{<:Integer}) = line_bundle(Z, degrees)

function twist(F::DegeneracyLocusBundle, i::Integer, k::Integer=1)
  X = base_variety(variety(F))
  1 <= i <= picard_rank(X) || throw(ArgumentError("Picard index is out of range."))
  degrees = zeros(Int, picard_rank(X))
  degrees[i] = k
  twist(F, degrees)
end
twist(F::DegeneracyLocusBundle, degrees::Vector{<:Integer}) =
  tensor_product(F, line_bundle(variety(F), degrees))

for op in (:structure_sheaf, :zero_bundle, :cotangent_bundle, :tangent_bundle,
  :canonical_bundle, :anticanonical_bundle)
  @eval $op(Z::DegeneracyLocus) = DegeneracyLocusBundle(Z, $op(incidence_resolution(Z)))
end
for op in (:dual, :det)
  @eval $op(F::DegeneracyLocusBundle) =
    DegeneracyLocusBundle(variety(F), $op(F.presentation))
end
determinant(F::DegeneracyLocusBundle) = det(F)
for op in (:exterior_power, :symmetric_power)
  @eval $op(F::DegeneracyLocusBundle, k::Integer) =
    DegeneracyLocusBundle(variety(F), $op(F.presentation, k))
end
for op in (:direct_sum, :tensor_product)
  @eval function $op(F::DegeneracyLocusBundle, H::DegeneracyLocusBundle)
    variety(F) == variety(H) || throw(ArgumentError("Bundles must live on the same locus."))
    DegeneracyLocusBundle(variety(F), $op(F.presentation, H.presentation))
  end
end
Base.:+(F::DegeneracyLocusBundle, H::DegeneracyLocusBundle) = direct_sum(F, H)
Base.:*(F::DegeneracyLocusBundle, H::DegeneracyLocusBundle) = tensor_product(F, H)
tensor_product(
  F::DegeneracyLocusBundle, H::Union{_AmbientBundle,RelativeGrassmannianBundle}
) =
  tensor_product(F, restrict(variety(F), H))
tensor_product(
  H::Union{_AmbientBundle,RelativeGrassmannianBundle}, F::DegeneracyLocusBundle
) =
  tensor_product(restrict(variety(F), H), F)
function Base.:*(n::Integer, F::DegeneracyLocusBundle)
  n >= 0 || throw(ArgumentError("Multiplicity must be nonnegative."))
  result = zero_bundle(variety(F))
  for _ in 1:n
    result = direct_sum(result, F)
  end
  result
end

cohomology(F::DegeneracyLocusBundle; kwargs...) = cohomology(F.presentation; kwargs...)
cohomology(Z::DegeneracyLocus, F::Union{_AmbientBundle,RelativeGrassmannianBundle}) =
  cohomology(restrict(Z, F))
chi(F::DegeneracyLocusBundle) = chi(F.presentation)
euler_characteristic(F::DegeneracyLocusBundle) = euler_characteristic(F.presentation)
degree(F::DegeneracyLocusBundle) = degree(F.presentation)
hodge_numbers(Z::DegeneracyLocus; kwargs...) =
  hodge_numbers(incidence_resolution(Z); kwargs...)
hochschild_cohomology(Z::DegeneracyLocus; kwargs...) =
  hochschild_cohomology(incidence_resolution(Z); kwargs...)
euler_characteristic(Z::DegeneracyLocus) = euler_characteristic(incidence_resolution(Z))

"""
    is_strongly_fano(Z::DegeneracyLocus) -> Bool

A sufficient anticanonical positivity certificate for the degeneracy locus.
The expected deeper stratum must be absent. The query then checks for an
ample lift on the incidence ambient, or an ample line on the base descending
from the incidence anticanonical bundle.

As for `ZeroLocus`, smoothness, nonemptiness and regularity are assumptions
about the section. A false result means no certificate was found; it does
not prove that the locus is not Fano.
"""
function is_strongly_fano(Z::DegeneracyLocus)
  d = dimension(base_variety(Z)) - rank(Z.data.base_equations)
  Z.data.lower_codimension > d || return false
  R = incidence_resolution(Z)
  is_strongly_fano(R) && return true
  if R isa ZeroLocus && Z.construction == :determinantal
    R = determinantal_resolution(Z.data.E, Z.data.F, Z.data.r;
      presentation=:relative, base_equations=Z.data.base_equations)
  end
  R isa RelativeZeroLocus || return false
  G = ambient_variety(R)
  L = tensor_product(anticanonical_bundle(G), dual(det(defining_bundle(R))))
  images = pushforward(G, L)
  rank(images[0]) == 1 && is_ample_line_bundle(_relative_total(images[0]))
end
