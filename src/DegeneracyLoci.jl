export BundleResolution, eagon_northcott, buchsbaum_rim, pfaffian_resolution
export gulliksen_negard

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
