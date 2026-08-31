# ═══════════════════════════════════════════════════════════════════════════════
#  Bundles on zero loci
#
#  A ZeroLocusBundle is represented by the terms of a bounded ambient
#  presentation. Its implicit differential makes the presentation exact away
#  from degree zero, where its cohomology is the represented bundle. These
#  presentations are closed under direct sums, tensor products and duals. In
#  characteristic zero, derived exterior and symmetric powers are computed in
#  the graded symmetric monoidal category of complexes.
# ═══════════════════════════════════════════════════════════════════════════════

export restrict

const _AmbientBundle = Union{CompletelyReducibleBundle,FilteredBundle}

"""Terms of a bounded ambient presentation which is exact away from degree zero."""
struct _AmbientBundlePresentation
  terms::Dict{Int,Vector{_AmbientBundle}}

  function _AmbientBundlePresentation(terms::Dict{Int,Vector{_AmbientBundle}})
    cleaned = Dict{Int,Vector{_AmbientBundle}}()
    for (degree, summands) in terms
      kept = _AmbientBundle[summand for summand in summands if rank(summand) != 0]
      isempty(kept) || (cleaned[degree] = kept)
    end
    new(cleaned)
  end
end

_AmbientBundlePresentation(degree::Int, bundle::_AmbientBundle) =
  _AmbientBundlePresentation(
    Dict{Int,Vector{_AmbientBundle}}(degree => _AmbientBundle[bundle])
  )

function _same_summands(left::Vector{_AmbientBundle}, right::Vector{_AmbientBundle})
  length(left) == length(right) || return false
  unmatched = trues(length(right))
  for summand in left
    index = findfirst(i -> unmatched[i] && summand == right[i], eachindex(right))
    index === nothing && return false
    unmatched[index] = false
  end
  true
end

function Base.:(==)(
  left::_AmbientBundlePresentation, right::_AmbientBundlePresentation
)
  keys(left.terms) == keys(right.terms) || return false
  all(
    degree -> _same_summands(left.terms[degree], right.terms[degree]),
    keys(left.terms),
  )
end

Base.isequal(
  left::_AmbientBundlePresentation, right::_AmbientBundlePresentation
) = left == right

function Base.hash(presentation::_AmbientBundlePresentation, h::UInt)
  term_hash = sum(
    hash((degree, sum(hash, summands; init=zero(UInt)))) for
    (degree, summands) in presentation.terms;
    init=zero(UInt),
  )
  hash(term_hash, h)
end

"""
    ZeroLocusBundle

A vector bundle on a [`ZeroLocus`](@ref), represented by the terms of a bounded
ambient presentation which is exact away from degree zero. The presentation
maps are implicit: the package records the data needed for bundle operations,
additive invariants, and symbolic long-exact-sequence calculations.

Use [`restrict`](@ref) for restrictions of ambient bundles and constructors
such as [`tangent_bundle`](@ref), [`cotangent_bundle`](@ref),
[`normal_bundle`](@ref), and [`conormal_bundle`](@ref) for intrinsic bundles.
"""
struct ZeroLocusBundle <: Bundle
  locus::ZeroLocus
  presentation::_AmbientBundlePresentation
end

Base.:(==)(F::ZeroLocusBundle, G::ZeroLocusBundle) =
  variety(F) == variety(G) && F.presentation == G.presentation
Base.isequal(F::ZeroLocusBundle, G::ZeroLocusBundle) = F == G
Base.hash(F::ZeroLocusBundle, h::UInt) = hash(F.presentation, hash(F.locus, h))

"""Return the zero locus on which `F` lives."""
variety(F::ZeroLocusBundle) = F.locus

"""Return the ambient partial flag variety of the base of `F`."""
ambient_variety(F::ZeroLocusBundle) = ambient_variety(variety(F))

function _check_same_locus(F::ZeroLocusBundle, G::ZeroLocusBundle, operation::String)
  variety(F) == variety(G) || throw(
    ArgumentError("$operation requires bundles on the same zero locus.")
  )
  nothing
end

"""
    restrict(Z::ZeroLocus, F) -> ZeroLocusBundle

Restrict an ambient [`CompletelyReducibleBundle`](@ref) or
[`FilteredBundle`](@ref) `F` to `Z`.

If `F` already lives on a formal zero locus cut out by `E`, it can be restricted
again to a zero locus in the same ambient variety whose defining bundle
contains `E` as a direct summand. Thus a bundle on `Z(E)` restricts naturally to
`Z(E ⊕ E′)`. Since [`ZeroLocus`](@ref) does not store individual sections,
containment is understood at the level of defining bundles.
"""
function restrict(Z::ZeroLocus, F::_AmbientBundle)
  marked_dynkin_type(variety(F)) == marked_dynkin_type(ambient_variety(Z)) || throw(
    ArgumentError("restrict requires a bundle on the ambient variety of the zero locus.")
  )
  ZeroLocusBundle(Z, _AmbientBundlePresentation(0, F))
end

function restrict(Z::ZeroLocus, F::ZeroLocusBundle)
  source = variety(F)
  source === Z && return F
  marked_dynkin_type(ambient_variety(source)) ==
  marked_dynkin_type(ambient_variety(Z)) || throw(
    ArgumentError("restrict requires zero loci in the same ambient variety type.")
  )
  is_summand(defining_bundle(source), defining_bundle(Z)) || throw(
    ArgumentError(
      "The target zero locus must contain the source defining bundle as a direct summand."
    ),
  )
  ZeroLocusBundle(Z, _AmbientBundlePresentation(F.presentation.terms))
end

"""Return the rank of a bundle represented by an ambient presentation."""
function rank(F::ZeroLocusBundle)
  result = sum(
    (
      (isodd(degree) ? -1 : 1) * sum(rank, summands; init=0) for
      (degree, summands) in F.presentation.terms
    );
    init=0,
  )
  result >= 0 || error("The ambient presentation has negative virtual rank $result.")
  result
end

function direct_sum(F::ZeroLocusBundle, G::ZeroLocusBundle)
  _check_same_locus(F, G, "direct_sum")
  terms = Dict{Int,Vector{_AmbientBundle}}(
    degree => copy(summands) for (degree, summands) in F.presentation.terms
  )
  for (degree, summands) in G.presentation.terms
    append!(get!(terms, degree, _AmbientBundle[]), summands)
  end
  ZeroLocusBundle(variety(F), _AmbientBundlePresentation(terms))
end

Base.:+(F::ZeroLocusBundle, G::ZeroLocusBundle) = direct_sum(F, G)

function tensor_product(F::ZeroLocusBundle, G::ZeroLocusBundle)
  _check_same_locus(F, G, "tensor_product")
  terms = Dict{Int,Vector{_AmbientBundle}}()
  for (degree_F, summands_F) in F.presentation.terms
    for (degree_G, summands_G) in G.presentation.terms
      degree = degree_F + degree_G
      degree_summands = get!(terms, degree, _AmbientBundle[])
      for summand_F in summands_F, summand_G in summands_G
        push!(degree_summands, tensor_product(summand_F, summand_G))
      end
    end
  end
  ZeroLocusBundle(variety(F), _AmbientBundlePresentation(terms))
end

Base.:*(F::ZeroLocusBundle, G::ZeroLocusBundle) = tensor_product(F, G)

function dual(F::ZeroLocusBundle)
  terms = Dict{Int,Vector{_AmbientBundle}}(
    -degree => _AmbientBundle[dual(summand) for summand in summands] for
    (degree, summands) in F.presentation.terms
  )
  ZeroLocusBundle(variety(F), _AmbientBundlePresentation(terms))
end

"""
    _derived_power(F, k, even_power, odd_power) -> ZeroLocusBundle

Construct a derived power of the ambient presentation of `F`. Each ambient
summand is paired with its cohomological degree. A multiexponent distributes
the total exponent `k` among those graded summands, and the resulting factors
are tensored in the sum of their weighted degrees.

Koszul signs exchange exterior and symmetric powers in odd degree. For a
derived exterior power, `even_power` is therefore `exterior_power` and
`odd_power` is `symmetric_power`; the derived symmetric power uses the reverse
assignment. The induced differentials are implicit, just as in the original
ambient presentation.
"""
function _derived_power(F::ZeroLocusBundle, k::Int, even_power, odd_power)
  k == 0 && return structure_sheaf(variety(F))

  presentation = F.presentation
  graded_summands = Tuple{Int,_AmbientBundle}[
    (degree, summand) for degree in sort!(collect(keys(presentation.terms))) for
    summand in presentation.terms[degree]
  ]
  isempty(graded_summands) && return zero_bundle(variety(F))

  terms = Dict{Int,Vector{_AmbientBundle}}()
  for multiplicities in multiexponents(length(graded_summands), k)
    factors = _AmbientBundle[]
    target_degree = 0
    for (multiplicity, (degree, summand)) in zip(multiplicities, graded_summands)
      multiplicity == 0 && continue
      target_degree += multiplicity * degree
      power_operation = iseven(degree) ? even_power : odd_power
      factor = power_operation(summand, multiplicity)
      rank(factor) == 0 && (empty!(factors); break)
      push!(factors, factor)
    end
    isempty(factors) && continue
    term = reduce(tensor_product, factors)
    push!(get!(terms, target_degree, _AmbientBundle[]), term)
  end
  ZeroLocusBundle(variety(F), _AmbientBundlePresentation(terms))
end

"""
    exterior_power(F::ZeroLocusBundle, k::Integer) -> ZeroLocusBundle

The `k`-th exterior power. For an arbitrary ambient presentation this uses
the derived exterior-power presentation, with symmetric powers on odd-degree
terms and exterior powers on even-degree terms. When `k` is greater than half
the rank, the perfect-pairing identity computes the complementary exterior
power of the dual and tensors it with [`det`](@ref).
"""
function exterior_power(F::ZeroLocusBundle, k::Integer)
  k = Int(k)
  r = rank(F)
  (k < 0 || k > r) && return zero_bundle(variety(F))
  k == 1 && return F
  k == r && return det(F)
  2k > r && return tensor_product(det(F), exterior_power(dual(F), r - k))
  _derived_power(F, k, exterior_power, symmetric_power)
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
  _derived_power(F, k, symmetric_power, exterior_power)
end

"""
    det(F::ZeroLocusBundle) -> ZeroLocusBundle

Compute the determinant line directly from the ambient presentation. A
summand in even cohomological degree contributes its determinant, while a
summand in odd degree contributes the dual of its determinant.
"""
function det(F::ZeroLocusBundle)
  factors = CompletelyReducibleBundle[]
  for degree in sort!(collect(keys(F.presentation.terms)))
    summands = F.presentation.terms[degree]
    for summand in summands
      factor = det(summand)
      push!(factors, isodd(degree) ? dual(factor) : factor)
    end
  end
  ambient_determinant =
    isempty(factors) ?
    structure_sheaf(ambient_variety(F)) :
    reduce(tensor_product, factors)
  restrict(variety(F), ambient_determinant)
end

determinant(F::ZeroLocusBundle) = det(F)

function Base.:*(n::Integer, F::ZeroLocusBundle)
  n < 0 && throw(ArgumentError("Cannot multiply a bundle by a negative integer ($n)"))
  n == 0 && return zero_bundle(variety(F))
  reduce(direct_sum, Iterators.repeated(F, n))
end

Base.:*(F::ZeroLocusBundle, n::Integer) = n * F
Base.iszero(F::ZeroLocusBundle) = iszero(rank(F))

"""
    structure_sheaf(Z::ZeroLocus) -> ZeroLocusBundle

Return the structure sheaf ``\\mathcal{O}_Z``.
"""
structure_sheaf(Z::ZeroLocus) = restrict(Z, structure_sheaf(ambient_variety(Z)))
O(Z::ZeroLocus) = structure_sheaf(Z)

"""
    zero_bundle(Z::ZeroLocus) -> ZeroLocusBundle

Return the rank-zero vector bundle on `Z`.
"""
zero_bundle(Z::ZeroLocus) =
  ZeroLocusBundle(
    Z, _AmbientBundlePresentation(Dict{Int,Vector{_AmbientBundle}}())
  )

"""
    line_bundle(Z::ZeroLocus, degree::Integer) -> ZeroLocusBundle
    line_bundle(Z::ZeroLocus, degrees::Vector{<:Integer}) -> ZeroLocusBundle

Restrict the ambient line bundle with the specified Picard degree or degrees
to `Z`.
"""
line_bundle(Z::ZeroLocus, degree::Integer) =
  restrict(Z, line_bundle(ambient_variety(Z), degree))
line_bundle(Z::ZeroLocus, degrees::Vector{<:Integer}) =
  restrict(Z, line_bundle(ambient_variety(Z), degrees))
O(Z::ZeroLocus, degree::Integer) = line_bundle(Z, degree)
O(Z::ZeroLocus, degrees::Vector{<:Integer}) = line_bundle(Z, degrees)

"""
    normal_bundle(Z::ZeroLocus) -> ZeroLocusBundle

Return the normal bundle ``\\mathcal{N}_{Z/X} = \\mathcal{E}|_Z`` of `Z` in its
ambient partial flag variety, where ``\\mathcal{E}`` is the defining bundle.
"""
normal_bundle(Z::ZeroLocus) = restrict(Z, defining_bundle(Z))

"""
    conormal_bundle(Z::ZeroLocus) -> ZeroLocusBundle

Return the conormal bundle ``\\mathcal{N}_{Z/X}^\\vee``.
"""
conormal_bundle(Z::ZeroLocus) = dual(normal_bundle(Z))

"""
    tangent_bundle(Z::ZeroLocus) -> ZeroLocusBundle

Return the tangent bundle of a smooth formal zero locus. Its ambient
presentation is the normal sequence
``0 \\to \\mathrm{T}_Z \\to \\mathrm{T}_X|_Z \\to \\mathcal{E}|_Z \\to 0``.
"""
function tangent_bundle(Z::ZeroLocus)
  terms = Dict{Int,Vector{_AmbientBundle}}(
    0 => _AmbientBundle[filtered_tangent_bundle(ambient_variety(Z))],
    1 => _AmbientBundle[defining_bundle(Z)],
  )
  ZeroLocusBundle(Z, _AmbientBundlePresentation(terms))
end

T(Z::ZeroLocus) = tangent_bundle(Z)

"""
    cotangent_bundle(Z::ZeroLocus) -> ZeroLocusBundle

Return the cotangent bundle of a smooth formal zero locus, represented by the
conormal sequence dual to the tangent presentation.
"""
function cotangent_bundle(Z::ZeroLocus)
  dual(tangent_bundle(Z))
end

"""
    canonical_bundle(Z::ZeroLocus) -> ZeroLocusBundle

Return the canonical bundle computed by adjunction,
``\\omega_Z = (\\omega_X \\otimes \\det \\mathcal{E})|_Z``.
"""
function canonical_bundle(Z::ZeroLocus)
  X = ambient_variety(Z)
  restrict(Z, tensor_product(canonical_bundle(X), det(defining_bundle(Z))))
end

"""Return the dual of [`canonical_bundle(Z)`](@ref)."""
anticanonical_bundle(Z::ZeroLocus) = dual(canonical_bundle(Z))

"""
    twist(F::ZeroLocusBundle, i::Integer, k::Integer=1) -> ZeroLocusBundle

Twist `F` by `O(k)` at the `i`-th marked node of its ambient partial flag
variety, matching the ambient-bundle interface.
"""
function twist(F::ZeroLocusBundle, i::Integer, k::Integer=1)
  X = ambient_variety(F)
  1 <= i <= picard_rank(X) || throw(
    ArgumentError(
      "Index $i out of range. The ambient has $(picard_rank(X)) marked node(s)."
    ),
  )
  degrees = zeros(Int, picard_rank(X))
  degrees[Int(i)] = Int(k)
  tensor_product(F, line_bundle(variety(F), degrees))
end

function twist(F::ZeroLocusBundle, degrees::Vector{<:Integer})
  tensor_product(F, line_bundle(variety(F), degrees))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Additive invariants and cohomology
# ═══════════════════════════════════════════════════════════════════════════════

_euler_characteristic_on_restriction(Z::ZeroLocus, F::CompletelyReducibleBundle) =
  _euler_characteristic_from_counts(Z, _to_counts(F))

function _euler_characteristic_on_restriction(Z::ZeroLocus, F::FilteredBundle)
  result = BigInt(0)
  for piece in graded_pieces(F)
    result += _euler_characteristic_on_restriction(Z, piece)
  end
  result
end

"""
    euler_characteristic(F::ZeroLocusBundle) -> BigInt

Compute the Euler characteristic of `F` as the alternating sum of the Euler
characteristics of its ambient presentation. This is exact even when the
individual cohomology groups are not determined by the induced maps.
"""
function euler_characteristic(F::ZeroLocusBundle)
  Z = variety(F)
  result = BigInt(0)
  for (degree, summands) in F.presentation.terms
    sign = isodd(degree) ? -1 : 1
    for summand in summands
      result += sign * _euler_characteristic_on_restriction(Z, summand)
    end
  end
  result
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
  Z::ZeroLocus, summands::Vector{_AmbientBundle}, var_counter::Ref{Int}
)
  isempty(summands) && return AffineExpr[AffineExpr(0) for _ in 0:dimension(Z)]
  entries = _restriction_cohomology(Z, first(summands), var_counter)
  for summand in Iterators.drop(summands, 1)
    summand_entries = _restriction_cohomology(Z, summand, var_counter)
    for i in eachindex(entries)
      entries[i] += summand_entries[i]
    end
  end
  entries
end

function _presentation_term_cohomology(
  Z::ZeroLocus, presentation::_AmbientBundlePresentation, degree::Int,
  var_counter::Ref{Int},
)
  _presentation_term_cohomology(
    Z, get(presentation.terms, degree, _AmbientBundle[]), var_counter
  )
end

"""
Compute cohomology from a complex which is exact away from degree zero.

The positive-degree half is reduced to its degree-zero kernel from right to
left. The negative-degree half is reduced to its image in degree zero from
left to right, and the requested bundle is the resulting cokernel.
"""
function _cohomology_from_presentation(F::ZeroLocusBundle, var_counter::Ref{Int})
  Z = variety(F)
  presentation = F.presentation
  isempty(presentation.terms) && return (
    AffineExpr[AffineExpr(0) for _ in 0:dimension(Z)],
    AffineExpr[],
  )

  minimum_degree = min(0, minimum(keys(presentation.terms)))
  maximum_degree = max(0, maximum(keys(presentation.terms)))
  inequalities = AffineExpr[]

  # K_0 = ker(C^0 -> C^1), computed through all positive degrees.
  kernel = _presentation_term_cohomology(
    Z, presentation, maximum_degree, var_counter
  )
  for degree in (maximum_degree - 1):-1:0
    term = _presentation_term_cohomology(Z, presentation, degree, var_counter)
    kernel = les_kernel(term, kernel, var_counter; inequalities)
  end

  minimum_degree == 0 && return (kernel, inequalities)

  # I_0 = im(C^-1 -> C^0), computed as the final cokernel of the negative
  # half. Then 0 -> I_0 -> K_0 -> F -> 0.
  image = _presentation_term_cohomology(Z, presentation, minimum_degree, var_counter)
  for degree in (minimum_degree + 1):-1
    term = _presentation_term_cohomology(Z, presentation, degree, var_counter)
    image = les_cokernel(image, term, var_counter; inequalities)
  end
  (les_cokernel(image, kernel, var_counter; inequalities), inequalities)
end

"""Return whether `F` is represented by the tangent normal sequence."""
function _is_tangent_bundle(F::ZeroLocusBundle)
  terms = F.presentation.terms
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

  # A degree-zero presentation is just a direct sum of restricted ambient
  # bundles, so no presentation LES or extra constraints are needed.
  if length(F.presentation.terms) == 1 && haskey(F.presentation.terms, 0)
    entries = _presentation_term_cohomology(Z, F.presentation.terms[0], var_counter)
    _renumber_variables!(entries)
    return Cohomology{AffineExpr}(entries, d)
  end

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

function Base.show(io::IO, F::ZeroLocusBundle)
  print(io, "Bundle(rank $(rank(F)) on ")
  show(io, variety(F))
  print(io, ")")
end
