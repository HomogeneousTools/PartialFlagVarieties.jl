export degeneracy_locus, determinantal_resolution, pfaffian_locus, symmetric_locus
export partially_decomposable_locus, rank_one_nilpotent_locus
export BundleResolution, eagon_northcott, buchsbaum_rim, pfaffian_resolution
export gulliksen_negard

function _check_degeneracy_data(E::_AmbientBundle, F::_AmbientBundle, r::Integer)
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
  E::_AmbientBundle, F::_AmbientBundle, r::Integer;
  presentation::Symbol=:auto, base_equations::_AmbientBundle=zero_bundle(variety(E)),
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
    degeneracy_locus(E, F, r; check_global_generation=false, check_smoothness=false, check_nonempty=false, kwargs...)

A rank-at-most-`r` degeneracy locus of a general map, retaining its incidence presentation.
The optional smoothness check certifies global generation and verifies that
expected codimension excludes the next rank stratum. Without that bound the
result is an incidence resolution whose invariants can differ from the image.
Failed sufficient certificates need not disprove smoothness.
No actual section is stored. Use `is_strongly_fano` on the result to query
the sufficient ambient anticanonical positivity certificate.
`check_nonempty=true` tests positive ample polarization degree under the
regularity assumptions; this optional intersection calculation can be costly.
"""
function degeneracy_locus(
  E::_AmbientBundle, F::_AmbientBundle, r::Integer;
  check_global_generation::Bool=false, check_smoothness::Bool=false,
  check_nonempty::Bool=false, kwargs...,
)
  _check_degeneracy_data(E, F, r)
  cuts = get(kwargs, :base_equations, zero_bundle(variety(E)))
  lower = r == 0 ? dimension(variety(E)) + 1 : (rank(E)-r+1)*(rank(F)-r+1)
  if check_global_generation || check_smoothness
    _check_degeneracy_geometry(tensor_product(dual(E), F), cuts, lower;
      check_global_generation, check_smoothness)
  end
  R = determinantal_resolution(E, F, r; check_nonempty, kwargs...)
  DegeneracyLocus(
    R, :determinantal, (; E, F, r=Int(r), base_equations=cuts, lower_codimension=lower)
  )
end

function _quadratic_locus(
  E::_AmbientBundle, L::CompletelyReducibleBundle, corank::Integer, symmetric::Bool;
  base_equations::_AmbientBundle=zero_bundle(variety(E)),
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
      check_global_generation, check_smoothness)
  end
  # At rank zero the form itself vanishes. The identity Grassmann bundle
  # also handles filtered equations without replacing their extensions.
  if corank == n
    G = GrassmannBundle(n, E)
    B = pullback(G, tensor_product(power(E, 2), L))
  else
    G = GrassmannBundle(n-Int(corank), E)
    S, Q = universal_subbundle(G), universal_quotient_bundle(G)
    pieces = [tensor_product(S, Q).representation, power(Q, 2).representation]
    W = RelativeGrassmannianBundle(G, FilteredBundle(G.auxiliary, pieces))
    B = tensor_product(W, pullback(G, L))
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
  E::_AmbientBundle, L::CompletelyReducibleBundle=structure_sheaf(variety(E));
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
  E::_AmbientBundle, L::CompletelyReducibleBundle=structure_sheaf(variety(E));
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
  E::_AmbientBundle; base_equations::_AmbientBundle=zero_bundle(variety(E)),
  check_global_generation::Bool=false, check_smoothness::Bool=false,
  check_nonempty::Bool=false,
)
  rank(E) == 6 || throw(ArgumentError("The three-form construction requires rank six."))
  variety(base_equations) == variety(E) || throw(ArgumentError("Incorrect base equations."))
  if check_global_generation || check_smoothness
    _check_degeneracy_geometry(exterior_power(E, 3), base_equations, 10;
      check_global_generation, check_smoothness)
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
  E::_AmbientBundle, L::CompletelyReducibleBundle;
  base_equations::_AmbientBundle=zero_bundle(variety(E)),
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
      check_global_generation, check_smoothness)
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
    tensor_product(tensor_product(dual(E), F), D),
    tensor_product(
      direct_sum(_trace_free_endomorphisms(E), _trace_free_endomorphisms(F)), D
    ),
    tensor_product(tensor_product(dual(F), E), D),
    tensor_product(D, D),
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
  for i in 1:(e - f + 1)
    push!(
      terms,
      tensor_product(
        tensor_product(exterior_power(E, f + i - 1), symmetric_power(dual(F), i - 1)), D
      ),
    )
  end
  BundleResolution(terms, dimension(variety(E)) - (e - f + 1))
end

"""Buchsbaum--Rim resolution of `coker(E → F)` for `rank(E) ≥ rank(F) > 0`, under the expected maximal-minor grade hypothesis."""
function buchsbaum_rim(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  variety(E) == variety(F) ||
    throw(ArgumentError("Source and target must live on the same base."))
  e, f = Int(rank(E)), Int(rank(F))
  e >= f > 0 || throw(ArgumentError("Require rank(E) ≥ rank(F) > 0."))
  terms = [F, E]
  for i in 2:(e - f + 1)
    push!(
      terms,
      tensor_product(
        tensor_product(exterior_power(E, f + i - 1), symmetric_power(dual(F), i - 2)),
        dual(det(F)),
      ),
    )
  end
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
    tensor_product(E, tensor_product(D, _line_power(L, -m))),
    tensor_product(dual(E), tensor_product(D, _line_power(L, -m - 1))),
    tensor_product(tensor_product(D, D), _line_power(L, -n))]
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
