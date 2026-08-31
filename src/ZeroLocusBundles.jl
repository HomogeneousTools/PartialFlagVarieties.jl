# ═══════════════════════════════════════════════════════════════════════════════
#  Bundles on zero loci
#
#  A ZeroLocusBundle is represented by a bounded complex of restrictions of
#  ambient bundles. The complex has a single cohomology sheaf, in degree zero,
#  equal to the bundle being represented. This presentation is closed under
#  direct sums, tensor products and duals. In characteristic zero, its derived
#  exterior and symmetric powers are computed in the graded symmetric monoidal
#  category of complexes.
# ═══════════════════════════════════════════════════════════════════════════════

export ZeroLocusBundle
export restrict, locus, presentation_degrees

"""A bounded complex of ambient bundles, exact away from degree zero."""
struct _AmbientBundleComplex
  terms::Dict{Int,Vector{Bundle}}

  function _AmbientBundleComplex(terms::Dict{Int,Vector{Bundle}})
    cleaned = Dict{Int,Vector{Bundle}}()
    for (degree, summands) in terms
      kept = Bundle[summand for summand in summands if rank(summand) != 0]
      isempty(kept) || (cleaned[degree] = kept)
    end
    new(cleaned)
  end
end

_AmbientBundleComplex(degree::Int, bundle::Bundle) =
  _AmbientBundleComplex(Dict{Int,Vector{Bundle}}(degree => Bundle[bundle]))

"""
    ZeroLocusBundle

A vector bundle on a [`ZeroLocus`](@ref), represented by a bounded complex of
restrictions of ambient bundles whose only cohomology sheaf is the represented
bundle in degree zero.

Use [`restrict`](@ref) for restrictions of ambient bundles and constructors
such as [`tangent_bundle`](@ref), [`cotangent_bundle`](@ref),
[`normal_bundle`](@ref), and [`conormal_bundle`](@ref) for intrinsic bundles.
"""
struct ZeroLocusBundle <: Bundle
  locus::ZeroLocus
  complex::_AmbientBundleComplex
  kind::Symbol
end

ZeroLocusBundle(Z::ZeroLocus, complex::_AmbientBundleComplex) =
  ZeroLocusBundle(Z, complex, :composite)

"""Return the zero locus on which `F` lives."""
locus(F::ZeroLocusBundle) = F.locus

"""Return the zero locus on which `F` lives."""
variety(F::ZeroLocusBundle) = locus(F)

"""Return the ambient partial flag variety of the base of `F`."""
ambient_variety(F::ZeroLocusBundle) = ambient_variety(locus(F))

"""Return the nonzero cohomological degrees in the ambient presentation."""
presentation_degrees(F::ZeroLocusBundle) = sort!(collect(keys(F.complex.terms)))

function _check_same_locus(F::ZeroLocusBundle, G::ZeroLocusBundle, operation::String)
  locus(F) === locus(G) || throw(
    ArgumentError("$operation requires bundles on the same zero locus.")
  )
  nothing
end

"""
    restrict(Z::ZeroLocus, F::Bundle) -> ZeroLocusBundle

Restrict an ambient bundle `F` to `Z`.
"""
function restrict(Z::ZeroLocus, F::Bundle)
  F isa ZeroLocusBundle && throw(
    ArgumentError("The bundle is already on a zero locus.")
  )
  marked_dynkin_type(variety(F)) == marked_dynkin_type(ambient_variety(Z)) || throw(
    ArgumentError("restrict requires a bundle on the ambient variety of the zero locus.")
  )
  ZeroLocusBundle(Z, _AmbientBundleComplex(0, F), :restriction)
end

"""Return the rank of a bundle represented by an ambient complex."""
function rank(F::ZeroLocusBundle)
  result = 0
  for (degree, summands) in F.complex.terms
    sign = isodd(degree) ? -1 : 1
    result += sign * sum(rank, summands; init=0)
  end
  result >= 0 || error("The ambient presentation has negative virtual rank $result.")
  result
end

function _direct_sum_complex(C::_AmbientBundleComplex, D::_AmbientBundleComplex)
  terms = Dict{Int,Vector{Bundle}}(
    degree => copy(summands) for (degree, summands) in C.terms
  )
  for (degree, summands) in D.terms
    append!(get!(terms, degree, Bundle[]), summands)
  end
  _AmbientBundleComplex(terms)
end

function direct_sum(F::ZeroLocusBundle, G::ZeroLocusBundle)
  _check_same_locus(F, G, "direct_sum")
  ZeroLocusBundle(locus(F), _direct_sum_complex(F.complex, G.complex))
end

Base.:+(F::ZeroLocusBundle, G::ZeroLocusBundle) = direct_sum(F, G)

function _tensor_complex(C::_AmbientBundleComplex, D::_AmbientBundleComplex)
  terms = Dict{Int,Vector{Bundle}}()
  for (degree_C, summands_C) in C.terms
    for (degree_D, summands_D) in D.terms
      degree = degree_C + degree_D
      target = get!(terms, degree, Bundle[])
      for F in summands_C, G in summands_D
        push!(target, tensor_product(F, G))
      end
    end
  end
  _AmbientBundleComplex(terms)
end

function tensor_product(F::ZeroLocusBundle, G::ZeroLocusBundle)
  _check_same_locus(F, G, "tensor_product")
  ZeroLocusBundle(locus(F), _tensor_complex(F.complex, G.complex))
end

Base.:*(F::ZeroLocusBundle, G::ZeroLocusBundle) = tensor_product(F, G)

function dual(F::ZeroLocusBundle)
  terms = Dict{Int,Vector{Bundle}}(
    -degree => Bundle[dual(summand) for summand in summands] for
    (degree, summands) in F.complex.terms
  )
  ZeroLocusBundle(locus(F), _AmbientBundleComplex(terms))
end

function _zero_locus_structure_sheaf(Z::ZeroLocus)
  restrict(Z, structure_sheaf(ambient_variety(Z)))
end

function _power_complex(C::_AmbientBundleComplex, k::Int, exterior::Bool, X)
  k < 0 && return _AmbientBundleComplex(Dict{Int,Vector{Bundle}}())
  k == 0 && return _AmbientBundleComplex(0, structure_sheaf(X))

  graded_summands = Tuple{Int,Bundle}[
    (degree, summand) for degree in sort!(collect(keys(C.terms))) for
    summand in C.terms[degree]
  ]
  isempty(graded_summands) && return _AmbientBundleComplex(Dict{Int,Vector{Bundle}}())

  terms = Dict{Int,Vector{Bundle}}()
  for alpha in multiexponents(length(graded_summands), k)
    factors = Bundle[]
    target_degree = 0
    for (multiplicity, (degree, summand)) in zip(alpha, graded_summands)
      multiplicity == 0 && continue
      target_degree += multiplicity * degree
      use_exterior = xor(exterior, isodd(degree))
      factor = if use_exterior
        exterior_power(summand, multiplicity)
      else
        symmetric_power(summand, multiplicity)
      end
      rank(factor) == 0 && (empty!(factors); break)
      push!(factors, factor)
    end
    isempty(factors) && continue
    term = reduce(tensor_product, factors)
    push!(get!(terms, target_degree, Bundle[]), term)
  end
  _AmbientBundleComplex(terms)
end

"""
    exterior_power(F::ZeroLocusBundle, k::Integer) -> ZeroLocusBundle

The `k`-th exterior power. For an arbitrary ambient presentation this uses
the derived exterior-power complex, with symmetric powers on odd-degree terms
and exterior powers on even-degree terms.
"""
function exterior_power(F::ZeroLocusBundle, k::Integer)
  k = Int(k)
  (k < 0 || k > rank(F)) && return zero_bundle(locus(F))
  k == 1 && return F
  ZeroLocusBundle(
    locus(F), _power_complex(F.complex, k, true, ambient_variety(F)), :exterior_power
  )
end

"""
    symmetric_power(F::ZeroLocusBundle, k::Integer) -> ZeroLocusBundle

The `k`-th symmetric power. For an arbitrary ambient presentation this uses
the derived symmetric-power complex, with exterior powers on odd-degree terms
and symmetric powers on even-degree terms.
"""
function symmetric_power(F::ZeroLocusBundle, k::Integer)
  k = Int(k)
  k < 0 && return zero_bundle(locus(F))
  k == 1 && return F
  ZeroLocusBundle(
    locus(F), _power_complex(F.complex, k, false, ambient_variety(F)), :symmetric_power
  )
end

det(F::ZeroLocusBundle) = exterior_power(F, rank(F))
determinant(F::ZeroLocusBundle) = det(F)

function Base.:*(n::Integer, F::ZeroLocusBundle)
  n < 0 && throw(ArgumentError("Cannot multiply a bundle by a negative integer ($n)"))
  n == 0 && return zero_bundle(locus(F))
  reduce(direct_sum, Iterators.repeated(F, n))
end

Base.:*(F::ZeroLocusBundle, n::Integer) = n * F
Base.iszero(F::ZeroLocusBundle) = rank(F) == 0 && isempty(F.complex.terms)

"""The structure sheaf of a zero locus."""
structure_sheaf(Z::ZeroLocus) = _zero_locus_structure_sheaf(Z)
O(Z::ZeroLocus) = structure_sheaf(Z)

"""The zero bundle on a zero locus."""
zero_bundle(Z::ZeroLocus) =
  ZeroLocusBundle(Z, _AmbientBundleComplex(Dict{Int,Vector{Bundle}}()), :zero)

line_bundle(Z::ZeroLocus, degree::Integer) =
  restrict(Z, line_bundle(ambient_variety(Z), degree))
line_bundle(Z::ZeroLocus, degrees::Vector{<:Integer}) =
  restrict(Z, line_bundle(ambient_variety(Z), degrees))
O(Z::ZeroLocus, degree::Integer) = line_bundle(Z, degree)
O(Z::ZeroLocus, degrees::Vector{<:Integer}) = line_bundle(Z, degrees)

"""The normal bundle of `Z` in its ambient partial flag variety."""
normal_bundle(Z::ZeroLocus) =
  ZeroLocusBundle(Z, _AmbientBundleComplex(0, defining_bundle(Z)), :normal)

"""The conormal bundle of `Z` in its ambient partial flag variety."""
conormal_bundle(Z::ZeroLocus) =
  ZeroLocusBundle(Z, _AmbientBundleComplex(0, dual(defining_bundle(Z))), :conormal)

"""The tangent bundle of a smooth zero locus, represented by its normal sequence."""
function tangent_bundle(Z::ZeroLocus)
  terms = Dict{Int,Vector{Bundle}}(
    0 => Bundle[filtered_tangent_bundle(ambient_variety(Z))],
    1 => Bundle[defining_bundle(Z)],
  )
  ZeroLocusBundle(Z, _AmbientBundleComplex(terms), :tangent)
end

T(Z::ZeroLocus) = tangent_bundle(Z)

"""The cotangent bundle of a smooth zero locus, represented by its conormal sequence."""
function cotangent_bundle(Z::ZeroLocus)
  F = dual(tangent_bundle(Z))
  ZeroLocusBundle(Z, F.complex, :cotangent)
end

"""The canonical bundle of a smooth zero locus, computed by adjunction."""
function canonical_bundle(Z::ZeroLocus)
  X = ambient_variety(Z)
  restrict(Z, tensor_product(canonical_bundle(X), det(defining_bundle(Z))))
end

"""The anticanonical bundle of a smooth zero locus."""
anticanonical_bundle(Z::ZeroLocus) = dual(canonical_bundle(Z))

function twist(F::ZeroLocusBundle, degree::Integer)
  tensor_product(F, line_bundle(locus(F), degree))
end

function twist(F::ZeroLocusBundle, degrees::Vector{<:Integer})
  tensor_product(F, line_bundle(locus(F), degrees))
end

function Base.show(io::IO, F::ZeroLocusBundle)
  print(io, "ZeroLocusBundle($(F.kind), rank $(rank(F)))")
end
