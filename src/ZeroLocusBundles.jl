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

export restrict

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
end

"""Return the zero locus on which `F` lives."""
variety(F::ZeroLocusBundle) = F.locus

"""Return the ambient partial flag variety of the base of `F`."""
ambient_variety(F::ZeroLocusBundle) = ambient_variety(variety(F))

function _check_same_locus(F::ZeroLocusBundle, G::ZeroLocusBundle, operation::String)
  variety(F) === variety(G) || throw(
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
  ZeroLocusBundle(Z, _AmbientBundleComplex(0, F))
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
  ZeroLocusBundle(variety(F), _direct_sum_complex(F.complex, G.complex))
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
  ZeroLocusBundle(variety(F), _tensor_complex(F.complex, G.complex))
end

Base.:*(F::ZeroLocusBundle, G::ZeroLocusBundle) = tensor_product(F, G)

function dual(F::ZeroLocusBundle)
  terms = Dict{Int,Vector{Bundle}}(
    -degree => Bundle[dual(summand) for summand in summands] for
    (degree, summands) in F.complex.terms
  )
  ZeroLocusBundle(variety(F), _AmbientBundleComplex(terms))
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
  (k < 0 || k > rank(F)) && return zero_bundle(variety(F))
  k == 1 && return F
  ZeroLocusBundle(variety(F), _power_complex(F.complex, k, true, ambient_variety(F)))
end

"""
    symmetric_power(F::ZeroLocusBundle, k::Integer) -> ZeroLocusBundle

The `k`-th symmetric power. For an arbitrary ambient presentation this uses
the derived symmetric-power complex, with exterior powers on odd-degree terms
and symmetric powers on even-degree terms.
"""
function symmetric_power(F::ZeroLocusBundle, k::Integer)
  k = Int(k)
  k < 0 && return zero_bundle(variety(F))
  k == 1 && return F
  ZeroLocusBundle(variety(F), _power_complex(F.complex, k, false, ambient_variety(F)))
end

det(F::ZeroLocusBundle) = exterior_power(F, rank(F))
determinant(F::ZeroLocusBundle) = det(F)

function Base.:*(n::Integer, F::ZeroLocusBundle)
  n < 0 && throw(ArgumentError("Cannot multiply a bundle by a negative integer ($n)"))
  n == 0 && return zero_bundle(variety(F))
  reduce(direct_sum, Iterators.repeated(F, n))
end

Base.:*(F::ZeroLocusBundle, n::Integer) = n * F
Base.iszero(F::ZeroLocusBundle) = rank(F) == 0 && isempty(F.complex.terms)

"""The structure sheaf of a zero locus."""
structure_sheaf(Z::ZeroLocus) = _zero_locus_structure_sheaf(Z)
O(Z::ZeroLocus) = structure_sheaf(Z)

"""The zero bundle on a zero locus."""
zero_bundle(Z::ZeroLocus) =
  ZeroLocusBundle(Z, _AmbientBundleComplex(Dict{Int,Vector{Bundle}}()))

line_bundle(Z::ZeroLocus, degree::Integer) =
  restrict(Z, line_bundle(ambient_variety(Z), degree))
line_bundle(Z::ZeroLocus, degrees::Vector{<:Integer}) =
  restrict(Z, line_bundle(ambient_variety(Z), degrees))
O(Z::ZeroLocus, degree::Integer) = line_bundle(Z, degree)
O(Z::ZeroLocus, degrees::Vector{<:Integer}) = line_bundle(Z, degrees)

"""The normal bundle of `Z` in its ambient partial flag variety."""
normal_bundle(Z::ZeroLocus) =
  ZeroLocusBundle(Z, _AmbientBundleComplex(0, defining_bundle(Z)))

"""The conormal bundle of `Z` in its ambient partial flag variety."""
conormal_bundle(Z::ZeroLocus) =
  ZeroLocusBundle(Z, _AmbientBundleComplex(0, dual(defining_bundle(Z))))

"""The tangent bundle of a smooth zero locus, represented by its normal sequence."""
function tangent_bundle(Z::ZeroLocus)
  terms = Dict{Int,Vector{Bundle}}(
    0 => Bundle[filtered_tangent_bundle(ambient_variety(Z))],
    1 => Bundle[defining_bundle(Z)],
  )
  ZeroLocusBundle(Z, _AmbientBundleComplex(terms))
end

T(Z::ZeroLocus) = tangent_bundle(Z)

"""The cotangent bundle of a smooth zero locus, represented by its conormal sequence."""
function cotangent_bundle(Z::ZeroLocus)
  dual(tangent_bundle(Z))
end

"""The canonical bundle of a smooth zero locus, computed by adjunction."""
function canonical_bundle(Z::ZeroLocus)
  X = ambient_variety(Z)
  restrict(Z, tensor_product(canonical_bundle(X), det(defining_bundle(Z))))
end

"""The anticanonical bundle of a smooth zero locus."""
anticanonical_bundle(Z::ZeroLocus) = dual(canonical_bundle(Z))

function twist(F::ZeroLocusBundle, degree::Integer)
  tensor_product(F, line_bundle(variety(F), degree))
end

function twist(F::ZeroLocusBundle, degrees::Vector{<:Integer})
  tensor_product(F, line_bundle(variety(F), degrees))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Additive invariants and cohomology
# ═══════════════════════════════════════════════════════════════════════════════

_euler_characteristic_on_restriction(Z::ZeroLocus, F::CompletelyReducibleBundle) =
  euler_characteristic(Z, F)

function _euler_characteristic_on_restriction(Z::ZeroLocus, F::FilteredBundle)
  sum(piece -> euler_characteristic(Z, piece), graded_pieces(F); init=BigInt(0))
end

"""
    euler_characteristic(F::ZeroLocusBundle) -> BigInt

Compute the Euler characteristic of `F` as the alternating sum of the Euler
characteristics of its ambient presentation. This is exact even when the
individual cohomology groups are not determined by the induced maps.
"""
function euler_characteristic(F::ZeroLocusBundle)
  Z = variety(F)
  sum(F.complex.terms; init=BigInt(0)) do (degree, summands)
    sign = isodd(degree) ? -1 : 1
    sign * sum(
      summand -> _euler_characteristic_on_restriction(Z, summand),
      summands;
      init=BigInt(0),
    )
  end
end

chi(F::ZeroLocusBundle) = euler_characteristic(F)

function _restriction_cohomology(
  Z::ZeroLocus, F::CompletelyReducibleBundle, var_counter::Ref{Int}
)
  _restrict_to_zero_locus_les(Z, _to_counts(F), var_counter)
end

function _restriction_cohomology(
  Z::ZeroLocus, F::FilteredBundle, var_counter::Ref{Int}
)
  _restrict_to_zero_locus_les(Z, F, var_counter)
end

function _presentation_term_cohomology(
  Z::ZeroLocus, summands::Vector{Bundle}, var_counter::Ref{Int}
)
  entries = AffineExpr[AffineExpr(0) for _ in 0:dimension(Z)]
  for summand in summands
    entries .+= _restriction_cohomology(Z, summand, var_counter)
  end
  entries
end

function _presentation_term_cohomology(
  Z::ZeroLocus, C::_AmbientBundleComplex, degree::Int, var_counter::Ref{Int}
)
  _presentation_term_cohomology(Z, get(C.terms, degree, Bundle[]), var_counter)
end

"""
Compute cohomology from a complex which is exact away from degree zero.

The positive-degree half is reduced to its degree-zero kernel from right to
left. The negative-degree half is reduced to its image in degree zero from
left to right, and the requested bundle is the resulting cokernel.
"""
function _cohomology_from_presentation(F::ZeroLocusBundle, var_counter::Ref{Int})
  Z = variety(F)
  C = F.complex
  isempty(C.terms) && return (
    AffineExpr[AffineExpr(0) for _ in 0:dimension(Z)],
    AffineExpr[],
  )

  minimum_degree = min(0, minimum(keys(C.terms)))
  maximum_degree = max(0, maximum(keys(C.terms)))
  inequalities = AffineExpr[]

  # K_0 = ker(C^0 -> C^1), computed through all positive degrees.
  kernel = _presentation_term_cohomology(Z, C, maximum_degree, var_counter)
  for degree in (maximum_degree - 1):-1:0
    term = _presentation_term_cohomology(Z, C, degree, var_counter)
    kernel = les_kernel(term, kernel, var_counter; inequalities)
  end

  minimum_degree == 0 && return (kernel, inequalities)

  # I_0 = im(C^-1 -> C^0), computed as the final cokernel of the negative
  # half. Then 0 -> I_0 -> K_0 -> F -> 0.
  image = _presentation_term_cohomology(Z, C, minimum_degree, var_counter)
  for degree in (minimum_degree + 1):-1
    term = _presentation_term_cohomology(Z, C, degree, var_counter)
    image = les_cokernel(image, term, var_counter; inequalities)
  end
  (les_cokernel(image, kernel, var_counter; inequalities), inequalities)
end

"""Return whether `F` is represented by the tangent normal sequence."""
function _is_tangent_bundle(F::ZeroLocusBundle)
  terms = F.complex.terms
  length(terms) == 2 && haskey(terms, 0) && haskey(terms, 1) || return false
  length(terms[0]) == 1 || return false
  length(terms[1]) == 1 || return false
  terms[0][1] == filtered_tangent_bundle(ambient_variety(F)) &&
    terms[1][1] == defining_bundle(variety(F))
end

"""
    cohomology(F::ZeroLocusBundle) -> Cohomology{AffineExpr}

Compute dimension-valued sheaf cohomology from the ambient presentation.
Entries are exact integers where exactness and the available geometric
constraints determine them, and symbolic affine expressions otherwise.

Character-valued cohomology is not defined: a section cutting out a zero locus
is generally not invariant under the ambient group.
"""
function cohomology(F::ZeroLocusBundle)
  Z = variety(F)
  d = dimension(Z)
  is_tangent = _is_tangent_bundle(F)

  if is_tangent && n_factors(Z) >= 2
    tangent_row = AffineExpr[hochschild_cohomology(Z)[1, q] for q in 0:d]
    all(is_determined, tangent_row) && return Cohomology{AffineExpr}(tangent_row, d)
  end

  var_counter = Ref(0)
  entries, inequalities = _cohomology_from_presentation(F, var_counter)
  entry_count = length(entries)
  system = vcat(entries, inequalities)

  # The presentation determines the Euler characteristic in K-theory.
  chi_equation =
    _alternating_sum(@view system[1:entry_count]) - AffineExpr(euler_characteristic(F))
  is_zero_expr(chi_equation) || _apply_equation!(system, chi_equation)
  _propagate_intervals!(system)

  # For a strict Calabi--Yau, H^0(T_Z) = H^0(Omega_Z^{d-1}) = 0.
  if is_tangent && d >= 2 && !is_zero_expr(system[1]) && is_strict_calabi_yau(Z)
    _apply_equation!(system, system[1])
    _propagate_intervals!(system)
  end

  entries = system[1:entry_count]
  _renumber_variables!(entries)
  Cohomology{AffineExpr}(entries, d)
end

"""Return whether every entry of a cohomology object is determined."""
is_determined(H::Cohomology) = all(is_determined, H.entries)

function Base.show(io::IO, F::ZeroLocusBundle)
  print(io, "Bundle(rank $(rank(F)) on ")
  show(io, variety(F))
  print(io, ")")
end
