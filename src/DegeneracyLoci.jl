export DegeneracyLocus, DegeneracyLocusBundle, incidence_resolution
export degeneracy_locus, determinantal_resolution, pfaffian_locus, symmetric_locus
export partially_decomposable_locus, rank_one_nilpotent_locus
export BundleResolution, eagon_northcott, buchsbaum_rim, pfaffian_resolution
export gulliksen_negard

# ─── Locus and bundle types ──────────────────────────────────────────────

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
  print(io, "Degeneracy locus(", Z.construction, ", dim ", dimension(Z),
    " in ", ambient_variety(Z), ")")

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

function restrict(Z::DegeneracyLocus, F::Union{AmbientBundle,RelativeGrassmannianBundle})
  R = incidence_resolution(Z)
  if R isa ZeroLocus && F isa AmbientBundle && variety(F) == base_variety(Z) &&
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
  F::DegeneracyLocusBundle, H::Union{AmbientBundle,RelativeGrassmannianBundle}
) =
  tensor_product(F, restrict(variety(F), H))
tensor_product(
  H::Union{AmbientBundle,RelativeGrassmannianBundle}, F::DegeneracyLocusBundle
) =
  tensor_product(restrict(variety(F), H), F)
function Base.:*(n::Integer, F::DegeneracyLocusBundle)
  n >= 0 || throw(ArgumentError("Multiplicity must be nonnegative."))
  reduce(direct_sum, Iterators.repeated(F, n); init=zero_bundle(variety(F)))
end

cohomology(F::DegeneracyLocusBundle; kwargs...) = cohomology(F.presentation; kwargs...)
cohomology(Z::DegeneracyLocus, F::Union{AmbientBundle,RelativeGrassmannianBundle}) =
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

# ─── Geometry certificates ───────────────────────────────────────────────

# A dominant graded piece is globally generated and has no higher cohomology.
# Induction through the filtration therefore lifts its generating sections.
_generation_certificate(F::CompletelyReducibleBundle) = is_globally_generated(F)
_generation_certificate(F::FilteredBundle) = all(is_globally_generated, graded_pieces(F))

function _require_generation(F::AmbientBundle, name::AbstractString)
  _generation_certificate(F) && return nothing
  reason = if F isa FilteredBundle
    "Global generation could not be certified"
  else
    "The bundle is not globally generated"
  end
  throw(ArgumentError("$reason for $name."))
end

function _check_degeneracy_geometry(
  section::AmbientBundle, cuts::AmbientBundle, lower_codimension::Integer;
  check_smoothness::Bool=false,
)
  variety(cuts) == variety(section) || throw(ArgumentError("Incorrect base equations."))
  _require_generation(section, "the section bundle")
  _require_generation(cuts, "the base equations")
  if check_smoothness
    d = dimension(variety(section)) - rank(cuts)
    lower_codimension > d || throw(
      ArgumentError("A deeper stratum may occur; smoothness could not be certified.")
    )
  end
  nothing
end

function _relative_generation_certificate(F::RelativeGrassmannianBundle)
  G = variety(F)
  b, n = rank(dynkin_type(base_variety(G))), Int(rank(G.bundle))
  for piece in _relative_pieces(F.representation)
    # The central GL(n) character is unrestricted; the other fibre coordinates
    # must be dominant for the relative evaluation map to be surjective.
    for rep in components(piece)
      c = coefficients(p_dominant_weight(rep))
      all(>=(0), c[(b + 1):(b + n - 1)]) || return false
    end
    images = pushforward(G, RelativeGrassmannianBundle(G, piece))
    _generation_certificate(images[0]) || return false
  end
  true
end

function _relative_ample_certificate(L::RelativeGrassmannianBundle)
  rank(L) == 1 || return false
  G = variety(L)
  if G.k == rank(G.bundle)
    return is_ample_line_bundle(_relative_total(pushforward(G, L)[0]))
  end
  R = _relative_total(L.representation)
  c = coefficients(p_dominant_weight(only(components(R))))
  c[rank(dynkin_type(base_variety(G))) + G.k] > 0 || return false
  H = line_bundle(base_variety(G), ones(Int, picard_rank(base_variety(G))))
  # Relative very ampleness plus generation gives a closed embedding into
  # X × P^N. Tensoring with an ample line on X is then ample.
  _relative_generation_certificate(tensor_product(L, pullback(G, dual(H))))
end

"""
    is_strongly_fano(Z::RelativeZeroLocus)

A sufficient certificate that the adjunction anticanonical line has an ample
lift to the relative ambient. A false result means the certificate failed;
it does not prove that the locus is not Fano. Regularity and nonemptiness
remain assumptions about the section.
"""
function is_strongly_fano(Z::RelativeZeroLocus)
  G = ambient_variety(Z)
  L = tensor_product(anticanonical_bundle(G), dual(det(defining_bundle(Z))))
  _relative_ample_certificate(L)
end

function _positive_polarization(G::GrassmannBundle)
  X = base_variety(G)
  H = line_bundle(X, ones(Int, picard_rank(X)))
  G.k == rank(G.bundle) && return pullback(G, H)
  images = pushforward(G, O(G, 1))
  # Twist every graded direct-image summand into the dominant chamber.
  t = 0
  for piece in _relative_pieces(images[0])
    for rep in components(piece), i in marked_nodes(X)
      t = max(t, -coefficients(p_dominant_weight(rep))[i])
    end
  end
  tensor_product(O(G, 1), pullback(G, _line_power(H, t + 1)))
end

function _check_degeneracy_nonempty(
  Z::Union{ZeroLocus,RelativeZeroLocus}; check_nonempty::Bool=false
)
  check_nonempty || return Z
  A = ambient_variety(Z)
  L = if A isa GrassmannBundle
    _positive_polarization(A)
  else
    line_bundle(A, ones(Int, picard_rank(A)))
  end
  volume = degree(restrict(Z, L))
  volume > 0 && return Z
  message = if iszero(volume)
    "The expected regular locus is empty."
  else
    "A negative polarization degree contradicts the regularity assumptions."
  end
  throw(ArgumentError(message))
end

# ─── Locus constructors and bundle resolutions ───────────────────────────

function _check_degeneracy_data(E::AmbientBundle, F::AmbientBundle, r::Integer)
  variety(E) == variety(F) ||
    throw(ArgumentError("Source and target must live on the same base."))
  0 <= r < min(rank(E), rank(F)) ||
    throw(ArgumentError("Require 0 ≤ r < min(rank(E), rank(F))."))
end

function _repeated_line(E::CompletelyReducibleBundle)
  isempty(components(E)) && return nothing
  u = first(components(E))
  degree(u) == 1 && all(==(u), components(E)) || return nothing
  CompletelyReducibleBundle(variety(E), [u])
end

"""
    determinantal_resolution(E, F, r; presentation=:auto, base_equations=..., check_global_generation=false, check_nonempty=false)

Kernel incidence resolution of the rank-at-most-`r` locus of a general map
`E → F`, with equations `S∨ ⊗ π*F` on `Gr(rank(E)-r,E)`.
The image may be singular; the result computes the incidence variety.
Regularity and expected codimension are assumed.

`:auto` recognizes products when either bundle is copies of one line.
`:relative` forces the relative presentation; `:product` requires a product.
`base_equations` may be filtered. Generation and nonemptiness checks are opt-in.
"""
function determinantal_resolution(
  E::AmbientBundle, F::AmbientBundle, r::Integer;
  presentation::Symbol=:auto, base_equations::AmbientBundle=zero_bundle(variety(E)),
  check_global_generation::Bool=false, check_nonempty::Bool=false,
)
  _check_degeneracy_data(E, F, r)
  X = variety(E)
  variety(base_equations) == X ||
    throw(ArgumentError("Base equations must live on the base."))
  presentation in (:auto, :product, :relative) ||
    throw(ArgumentError("Unknown presentation: $presentation."))
  H = check_global_generation || r == 0 ? tensor_product(dual(E), F) : nothing
  if check_global_generation
    _require_generation(H, "Hom(E,F)")
    _require_generation(base_equations, "the base equations")
  end
  if r == 0 && H isa CompletelyReducibleBundle &&
    base_equations isa CompletelyReducibleBundle
    return _check_degeneracy_nonempty(
      zero_locus(direct_sum(H, base_equations)); check_nonempty
    )
  end
  if r > 0 && presentation != :relative && E isa CompletelyReducibleBundle &&
    F isa CompletelyReducibleBundle && base_equations isa CompletelyReducibleBundle
    L = _repeated_line(E)
    if L === nothing && _repeated_line(F) !== nothing
      return determinantal_resolution(
        dual(F), dual(E), r; presentation, base_equations, check_nonempty
      )
    elseif L !== nothing
      A = Gr(Int(rank(E)) - Int(r), Int(rank(E)))
      B = external_tensor_product(tensor_product(dual(L), F), dual(universal_subbundle(A)))
      cuts = external_tensor_product(base_equations, structure_sheaf(A))
      return _check_degeneracy_nonempty(zero_locus(direct_sum(B, cuts)); check_nonempty)
    end
  end
  presentation == :product &&
    throw(ArgumentError("No completely reducible product presentation is available."))
  G = GrassmannBundle(Int(rank(E)) - Int(r), E)
  B = tensor_product(dual(universal_subbundle(G)), pullback(G, F))
  _check_degeneracy_nonempty(
    zero_locus(direct_sum(B, pullback(G, base_equations))); check_nonempty
  )
end

"""
    degeneracy_locus(E, F, r; presentation=:auto, base_equations=...,
                     check_global_generation=false, check_smoothness=false, check_nonempty=false)

A rank-at-most-`r` degeneracy locus of a general map, retaining its incidence presentation.
The optional smoothness check certifies global generation and verifies that
expected codimension excludes the next rank stratum. Without that bound the
result is an incidence resolution whose invariants can differ from the image.
Failed sufficient certificates need not disprove smoothness.
No actual section is stored.
`check_nonempty=true` tests positive ample polarization degree under the
regularity assumptions; this optional intersection calculation can be costly.
"""
function degeneracy_locus(
  E::AmbientBundle, F::AmbientBundle, r::Integer;
  presentation::Symbol=:auto, base_equations::AmbientBundle=zero_bundle(variety(E)),
  check_global_generation::Bool=false, check_smoothness::Bool=false,
  check_nonempty::Bool=false,
)
  _check_degeneracy_data(E, F, r)
  lower = r == 0 ? dimension(variety(E)) + 1 : (rank(E)-r+1)*(rank(F)-r+1)
  if check_global_generation || check_smoothness
    _check_degeneracy_geometry(tensor_product(dual(E), F), base_equations, lower;
      check_smoothness)
  end
  R = determinantal_resolution(E, F, r; presentation, base_equations, check_nonempty)
  DegeneracyLocus(
    R, :determinantal, (; E, F, r=Int(r), base_equations, lower_codimension=lower)
  )
end

function _quadratic_locus(
  E::AmbientBundle, L::CompletelyReducibleBundle, corank::Integer, symmetric::Bool;
  base_equations::AmbientBundle=zero_bundle(variety(E)),
  check_global_generation::Bool=false, check_smoothness::Bool=false,
  check_nonempty::Bool=false,
)
  X, n = variety(E), Int(rank(E))
  variety(L) == X && rank(L) == 1 ||
    throw(ArgumentError("L must be a line bundle on the base."))
  variety(base_equations) == X || throw(ArgumentError("Incorrect base equations."))
  n > 0 && 0 <= corank <= n ||
    throw(ArgumentError("Require 0 ≤ corank ≤ rank(E) and positive rank."))
  symmetric || iseven(n-corank) ||
    throw(ArgumentError("The allowed skew rank must be even."))
  power = symmetric ? symmetric_power : exterior_power
  next = Int(corank) + (symmetric ? 1 : 2)
  lower = next > n ? dimension(X)+1 : (symmetric ? next*(next+1)÷2 : next*(next-1)÷2)
  if check_global_generation || check_smoothness
    _check_degeneracy_geometry(tensor_product(power(E, 2), L), base_equations, lower;
      check_smoothness)
  end
  # At rank zero the form itself vanishes. The identity Grassmann bundle
  # also handles filtered equations without replacing their extensions.
  G = GrassmannBundle(corank == n ? n : n-Int(corank), E)
  B = if corank == n
    pullback(G, tensor_product(power(E, 2), L))
  else
    S, Q = universal_subbundle(G), universal_quotient_bundle(G)
    pieces = [tensor_product(S, Q).representation, power(Q, 2).representation]
    W = RelativeGrassmannianBundle(G, FilteredBundle(G.auxiliary, pieces))
    tensor_product(W, pullback(G, L))
  end
  Z = zero_locus(direct_sum(B, pullback(G, base_equations)))
  _check_degeneracy_nonempty(Z; check_nonempty)
  DegeneracyLocus(Z, symmetric ? :symmetric : :pfaffian,
    (; E, L, corank=Int(corank), base_equations, lower_codimension=lower))
end

"""
    pfaffian_locus(E, L=structure_sheaf(variety(E)); corank=3, kwargs...)

A degeneracy locus for a section of `∧²E ⊗ L` of corank at least
`corank`. The allowed rank `rank(E)-corank` must be even. The expected
codimension is `corank*(corank-1)/2`. Equations retain the filtration of
`∧²π*E / ∧²S` on `Gr(rank(E)-corank,E)`.
Global-generation, smoothness and nonemptiness checks are optional, as for
[`degeneracy_locus`](@ref). Without a smoothness certificate the result may
resolve a singular image. Coranks zero and full rank are also supported.
"""
function pfaffian_locus(
  E::AmbientBundle, L::CompletelyReducibleBundle=structure_sheaf(variety(E));
  corank::Integer=3, kwargs...,
)
  _quadratic_locus(E, L, corank, false; kwargs...)
end

"""
    symmetric_locus(E, L=structure_sheaf(variety(E)); corank=1, kwargs...)

A degeneracy locus for a section of `Sym²E ⊗ L` of corank at
least `corank`, with expected codimension `corank*(corank+1)/2`.
The equations are the filtered quotient `Sym²π*E / Sym²S` on
`Gr(rank(E)-corank,E)`. Accepts the same optional geometry checks as
[`degeneracy_locus`](@ref).
"""
function symmetric_locus(
  E::AmbientBundle, L::CompletelyReducibleBundle=structure_sheaf(variety(E));
  corank::Integer=1, kwargs...,
)
  _quadratic_locus(E, L, corank, true; kwargs...)
end

"""
    partially_decomposable_locus(E; base_equations=..., kwargs...)

Rank-six three-form construction from arXiv:1704.01436, Section 3 and
Appendix A, with equations `∧³Q` on `Gr(1,E)`.
Accepts optional global-generation, smoothness and nonemptiness checks.
The smoothness certificate excludes the decomposable stratum by dimension.
"""
function partially_decomposable_locus(
  E::AmbientBundle; base_equations::AmbientBundle=zero_bundle(variety(E)),
  check_global_generation::Bool=false, check_smoothness::Bool=false,
  check_nonempty::Bool=false,
)
  rank(E) == 6 || throw(ArgumentError("The three-form construction requires rank six."))
  variety(base_equations) == variety(E) || throw(ArgumentError("Incorrect base equations."))
  if check_global_generation || check_smoothness
    _check_degeneracy_geometry(exterior_power(E, 3), base_equations, 10;
      check_smoothness)
  end
  G = GrassmannBundle(1, E)
  Z = zero_locus(
    direct_sum(exterior_power(universal_quotient_bundle(G), 3), pullback(G, base_equations))
  )
  _check_degeneracy_nonempty(Z; check_nonempty)
  DegeneracyLocus(Z, :three_form, (; E, base_equations, lower_codimension=10))
end

"""
    rank_one_nilpotent_locus(E, L; base_equations=..., kwargs...)

Incidence presentation of traceless rank-at-most-one endomorphisms in
`End(E) ⊗ L` on `Gr(1,E)`. The allowed subbundle is `Hom(Q,S) ⊗ L`;
the quotient has graded pieces `End(Q) ⊗ L` and `Hom(S,Q) ⊗ L`.
Accepts optional global-generation, smoothness and nonemptiness checks.
The smoothness certificate excludes the zero endomorphism stratum.
"""
function rank_one_nilpotent_locus(
  E::AmbientBundle, L::CompletelyReducibleBundle;
  base_equations::AmbientBundle=zero_bundle(variety(E)),
  check_global_generation::Bool=false, check_smoothness::Bool=false,
  check_nonempty::Bool=false,
)
  variety(L) == variety(E) && rank(L) == 1 ||
    throw(ArgumentError("L must be a line bundle on the base."))
  variety(base_equations) == variety(E) || throw(ArgumentError("Incorrect base equations."))
  rank(E) >= 2 || throw(ArgumentError("Require rank at least two."))
  if check_global_generation || check_smoothness
    adjoint = tensor_product(E, dual(E))
    # Generation of End(E) suffices for its trace-free direct summand.
    _check_degeneracy_geometry(tensor_product(adjoint, L), base_equations, rank(E)^2-1;
      check_smoothness)
  end
  G = GrassmannBundle(1, E)
  S, Q = universal_subbundle(G), universal_quotient_bundle(G)
  pieces = [
    tensor_product(Q, dual(Q)).representation, tensor_product(dual(S), Q).representation
  ]
  W = RelativeGrassmannianBundle(G, FilteredBundle(G.auxiliary, pieces))
  Z = zero_locus(direct_sum(tensor_product(W, pullback(G, L)), pullback(G, base_equations)))
  _check_degeneracy_nonempty(Z; check_nonempty)
  DegeneracyLocus(Z, :nilpotent, (; E, L, base_equations, lower_codimension=rank(E)^2-1))
end

"""
    BundleResolution

Terms of a locally free resolution `… → K₂ → K₁ → K₀ → M → 0` on a
homogeneous base. Stores terms and support dimension, not differentials.
The constructors assume the expected grade conditions for exactness.
`cohomology` keeps undetermined induced-map ranks symbolic.
"""
struct BundleResolution
  terms::Vector{CompletelyReducibleBundle}
  support_dimension::Int
  function BundleResolution(terms::Vector{CompletelyReducibleBundle}, d::Integer)
    isempty(terms) && throw(ArgumentError("A resolution needs a degree-zero term."))
    X = variety(first(terms))
    all(variety(A) == X for A in terms) ||
      throw(ArgumentError("All resolution terms must live on the same base."))
    0 <= d <= dimension(X) || throw(ArgumentError("Invalid expected support dimension."))
    new(terms, Int(d))
  end
end
Base.getindex(R::BundleResolution, i::Int) = R.terms[i + 1]
Base.firstindex(::BundleResolution) = 0
Base.lastindex(R::BundleResolution) = length(R.terms) - 1

function _trace_free_endomorphisms(E::CompletelyReducibleBundle)
  X = variety(E)
  counts = _virtual_bundle(tensor_product(E, dual(E)))
  unit = only(components(structure_sheaf(X)))
  counts[unit] -= 1
  _genuine_bundle(X, counts)
end

"""
    gulliksen_negard(E, F)

Structure-sheaf resolution of the submaximal-minor locus of `E → F`, for
equal ranks at least two and expected codimension four. If
`D = det(E) ⊗ det(F)∨`, the positive-degree terms are `Hom(E,F) ⊗ D`,
`(sl(E) ⊕ sl(F)) ⊗ D`, `Hom(F,E) ⊗ D`, and `D²`.
The expected-grade hypothesis is assumed; no section is stored.
"""
function gulliksen_negard(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  variety(E) == variety(F) ||
    throw(ArgumentError("Source and target must live on the same base."))
  rank(E) == rank(F) >= 2 || throw(ArgumentError("Require equal ranks at least two."))
  D = tensor_product(det(E), dual(det(F)))
  terms = [
    structure_sheaf(variety(E)),
    dual(E) ⊗ F ⊗ D,
    (_trace_free_endomorphisms(E) + _trace_free_endomorphisms(F)) ⊗ D,
    dual(F) ⊗ E ⊗ D,
    D ⊗ D,
  ]
  BundleResolution(terms, dimension(variety(E)) - 4)
end

"""Eagon--Northcott resolution of the structure sheaf of the maximal-minor locus of `E → F`, with `rank(E) ≥ rank(F)` and expected codimension."""
function eagon_northcott(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  variety(E) == variety(F) ||
    throw(ArgumentError("Source and target must live on the same base."))
  e, f = Int(rank(E)), Int(rank(F))
  e >= f > 0 || throw(ArgumentError("Require rank(E) ≥ rank(F) > 0."))
  D = dual(det(F))
  terms = [structure_sheaf(variety(E))]
  append!(
    terms, [
      exterior_power(E, f+i) ⊗ symmetric_power(dual(F), i) ⊗ D for i in 0:(e - f)
    ]
  )
  BundleResolution(terms, dimension(variety(E)) - (e - f + 1))
end

"""Buchsbaum--Rim resolution of `coker(E → F)` for `rank(E) ≥ rank(F) > 0`, under the expected maximal-minor grade hypothesis."""
function buchsbaum_rim(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  variety(E) == variety(F) ||
    throw(ArgumentError("Source and target must live on the same base."))
  e, f = Int(rank(E)), Int(rank(F))
  e >= f > 0 || throw(ArgumentError("Require rank(E) ≥ rank(F) > 0."))
  D = dual(det(F))
  terms = [F, E]
  append!(
    terms, [
      exterior_power(E, f+i) ⊗ symmetric_power(dual(F), i-1) ⊗ D for i in 1:(e - f)
    ]
  )
  BundleResolution(terms, dimension(variety(E)) - (e - f + 1))
end

"""Buchsbaum--Eisenbud structure-sheaf resolution for the expected codimension-three Pfaffian locus of `∧²E ⊗ L`, with odd rank at least three."""
function pfaffian_resolution(
  E::CompletelyReducibleBundle, L::CompletelyReducibleBundle=structure_sheaf(variety(E))
)
  variety(L) == variety(E) && rank(L) == 1 ||
    throw(ArgumentError("L must be a line bundle on the base."))
  n = Int(rank(E))
  isodd(n) && n >= 3 || throw(ArgumentError("Require odd rank at least three."))
  m = (n - 1) ÷ 2
  D = dual(det(E))
  terms = [structure_sheaf(variety(E)),
    E ⊗ D ⊗ _line_power(L, -m),
    dual(E) ⊗ D ⊗ _line_power(L, -m - 1),
    D ⊗ D ⊗ _line_power(L, -n)]
  BundleResolution(terms, dimension(variety(E)) - 3)
end

function tensor_product(R::BundleResolution, F::CompletelyReducibleBundle)
  BundleResolution([tensor_product(A, F) for A in R.terms], R.support_dimension)
end
euler_characteristic(R::BundleResolution) = sum(
  ((-1)^(i - 1) * euler_characteristic(A) for (i, A) in enumerate(R.terms));
  init=BigInt(0),
)
chi(R::BundleResolution) = euler_characteristic(R)
function cohomology(R::BundleResolution)
  counter = Ref(0)
  inequalities = AffineExpr[]
  rows = [AffineExpr.(cohomology(A).entries) for A in reverse(R.terms)]
  entries = long_exact_sequence_cokernel(rows, counter; inequalities)
  n = length(entries)
  system = vcat(entries, inequalities)
  for i in (R.support_dimension + 2):n
    is_zero_expr(system[i]) || _apply_equation!(system, system[i])
  end
  _propagate_intervals!(system)
  entries = system[1:(R.support_dimension + 1)]
  _renumber_variables!(entries)
  Cohomology{AffineExpr}(entries, R.support_dimension)
end
