export RelativeZeroLocus, RelativeZeroLocusBundle

"""
    zero_locus(F::RelativeGrassmannianBundle; check_global_generation=false)

A regular zero locus in a relative Grassmannian. Smoothness and regularity
are assumptions about the chosen section, which is not stored.
The optional global-generation check uses a sufficient relative certificate.
"""
mutable struct RelativeZeroLocus
  const defining_bundle::RelativeGrassmannianBundle
  wedges::Union{Nothing,Vector{RelativeGrassmannianBundle}}
end

function zero_locus(F::RelativeGrassmannianBundle; check_global_generation::Bool=false)
  rank(F) <= dimension(variety(F)) ||
    throw(ArgumentError("Bundle rank exceeds ambient dimension."))
  check_global_generation && !_relative_generation_certificate(F) &&
    throw(
      ArgumentError(
        "Global generation of the relative defining bundle could not be certified."
      ),
    )
  RelativeZeroLocus(F, nothing)
end
defining_bundle(Z::RelativeZeroLocus) = Z.defining_bundle
ambient_variety(Z::RelativeZeroLocus) = variety(defining_bundle(Z))
dimension(Z::RelativeZeroLocus) =
  dimension(ambient_variety(Z)) - Int(rank(defining_bundle(Z)))
codimension(Z::RelativeZeroLocus) = Int(rank(defining_bundle(Z)))
Base.:(==)(Z::RelativeZeroLocus, W::RelativeZeroLocus) =
  defining_bundle(Z) == defining_bundle(W)
Base.isequal(Z::RelativeZeroLocus, W::RelativeZeroLocus) = Z == W
Base.hash(Z::RelativeZeroLocus, h::UInt) = hash(defining_bundle(Z), h)
Base.show(io::IO, Z::RelativeZeroLocus) =
  print(io, "Zero locus(dim ", dimension(Z), " in ", ambient_variety(Z), ")")

function _relative_wedges!(Z::RelativeZeroLocus)
  if Z.wedges === nothing
    D = dual(defining_bundle(Z))
    Z.wedges = [exterior_power(D, i) for i in 0:codimension(Z)]
  end
  Z.wedges::Vector{RelativeGrassmannianBundle}
end

"""A bundle on a relative zero locus, with a bounded ambient presentation."""
struct RelativeZeroLocusBundle <: Bundle
  variety::RelativeZeroLocus
  terms::Dict{Int,RelativeGrassmannianBundle}

  function RelativeZeroLocusBundle(
    Z::RelativeZeroLocus, terms::Dict{Int,RelativeGrassmannianBundle}
  )
    all(variety(F) == ambient_variety(Z) for F in values(terms)) ||
      throw(ArgumentError("Presentation terms must live on the ambient Grassmann bundle."))
    new(Z, filter(p -> !iszero(last(p)), terms))
  end
end
variety(F::RelativeZeroLocusBundle) = F.variety
ambient_variety(F::RelativeZeroLocusBundle) = ambient_variety(variety(F))
rank(F::RelativeZeroLocusBundle) = sum(((-1)^i * rank(A) for (i, A) in F.terms); init=0)
# A rank-zero vector bundle is zero even when its exact ambient presentation
# has nonzero terms, as for the tangent bundle of a reduced point.
Base.iszero(F::RelativeZeroLocusBundle) = iszero(rank(F))
Base.:(==)(F::RelativeZeroLocusBundle, H::RelativeZeroLocusBundle) =
  variety(F) == variety(H) && F.terms == H.terms
Base.isequal(F::RelativeZeroLocusBundle, H::RelativeZeroLocusBundle) = F == H
Base.hash(F::RelativeZeroLocusBundle, h::UInt) = hash(F.terms, hash(variety(F), h))
Base.show(io::IO, F::RelativeZeroLocusBundle) =
  print(io, "Bundle(rank ", rank(F), " on ", variety(F), ")")
restrict(Z::RelativeZeroLocus, F::RelativeGrassmannianBundle) =
  RelativeZeroLocusBundle(Z, Dict(0 => F))
restrict(Z::RelativeZeroLocus, F::_AmbientBundle) =
  restrict(Z, pullback(ambient_variety(Z), F))
structure_sheaf(Z::RelativeZeroLocus) = restrict(Z, structure_sheaf(ambient_variety(Z)))
zero_bundle(Z::RelativeZeroLocus) =
  RelativeZeroLocusBundle(Z, Dict{Int,RelativeGrassmannianBundle}())

function _relative_add_term!(terms, i, F)
  terms[i] = haskey(terms, i) ? direct_sum(terms[i], F) : F
end

function direct_sum(F::RelativeZeroLocusBundle, H::RelativeZeroLocusBundle)
  variety(F) == variety(H) || throw(ArgumentError("Bundles must live on the same locus."))
  terms = copy(F.terms)
  for (i, A) in H.terms
    _relative_add_term!(terms, i, A)
  end
  RelativeZeroLocusBundle(variety(F), terms)
end
function tensor_product(F::RelativeZeroLocusBundle, H::RelativeZeroLocusBundle)
  variety(F) == variety(H) || throw(ArgumentError("Bundles must live on the same locus."))
  terms = Dict{Int,RelativeGrassmannianBundle}()
  for (i, A) in F.terms, (j, B) in H.terms
    _relative_add_term!(terms, i + j, tensor_product(A, B))
  end
  RelativeZeroLocusBundle(variety(F), terms)
end
Base.:+(F::RelativeZeroLocusBundle, H::RelativeZeroLocusBundle) = direct_sum(F, H)
Base.:*(F::RelativeZeroLocusBundle, H::RelativeZeroLocusBundle) = tensor_product(F, H)
dual(F::RelativeZeroLocusBundle) =
  RelativeZeroLocusBundle(variety(F), Dict(-i => dual(A) for (i, A) in F.terms))

function _relative_derived_power(F::RelativeZeroLocusBundle, k::Int, even_power, odd_power)
  k < 0 && return zero_bundle(variety(F))
  k == 0 && return structure_sheaf(variety(F))
  isempty(F.terms) && return zero_bundle(variety(F))
  summands = sort!(collect(F.terms); by=first)
  terms = Dict{Int,RelativeGrassmannianBundle}()
  for α in multiexponents(length(summands), k)
    term = structure_sheaf(ambient_variety(F))
    degree = 0
    for (m, (i, A)) in zip(α, summands)
      m == 0 && continue
      term = tensor_product(term, (iseven(i) ? even_power : odd_power)(A, m))
      degree += m * i
      iszero(term) && break
    end
    iszero(term) || _relative_add_term!(terms, degree, term)
  end
  RelativeZeroLocusBundle(variety(F), terms)
end
function exterior_power(F::RelativeZeroLocusBundle, k::Integer)
  r = Int(rank(F))
  (k < 0 || k > r) && return zero_bundle(variety(F))
  k == 0 && return structure_sheaf(variety(F))
  k == 1 && return F
  2k > r && return tensor_product(det(F), exterior_power(dual(F), r - k))
  _relative_derived_power(F, Int(k), exterior_power, symmetric_power)
end
symmetric_power(F::RelativeZeroLocusBundle, k::Integer) =
  _relative_derived_power(F, Int(k), symmetric_power, exterior_power)
function det(F::RelativeZeroLocusBundle)
  D = structure_sheaf(ambient_variety(F))
  for (i, A) in F.terms
    D = tensor_product(D, iseven(i) ? det(A) : dual(det(A)))
  end
  restrict(variety(F), D)
end
determinant(F::RelativeZeroLocusBundle) = det(F)

normal_bundle(Z::RelativeZeroLocus) = restrict(Z, defining_bundle(Z))
conormal_bundle(Z::RelativeZeroLocus) = dual(normal_bundle(Z))
cotangent_bundle(Z::RelativeZeroLocus) = RelativeZeroLocusBundle(
  Z, Dict(-1 => dual(defining_bundle(Z)), 0 => cotangent_bundle(ambient_variety(Z)))
)
tangent_bundle(Z::RelativeZeroLocus) = dual(cotangent_bundle(Z))
canonical_bundle(Z::RelativeZeroLocus) =
  restrict(Z, tensor_product(canonical_bundle(ambient_variety(Z)), det(defining_bundle(Z))))
anticanonical_bundle(Z::RelativeZeroLocus) = dual(canonical_bundle(Z))

function _relative_restriction_chi(Z::RelativeZeroLocus, F::RelativeGrassmannianBundle)
  sum(
    (
      (-1)^i * euler_characteristic(tensor_product(F, W)) for
      (i, W) in enumerate(_relative_wedges!(Z))
    );
    init=BigInt(0),
  ) * -1
end
euler_characteristic(F::RelativeZeroLocusBundle) = sum(
  ((-1)^i * _relative_restriction_chi(variety(F), A) for (i, A) in F.terms);
  init=BigInt(0),
)
chi(F::RelativeZeroLocusBundle) = euler_characteristic(F)

function _relative_restriction_cohomology(
  Z::RelativeZeroLocus, F::RelativeGrassmannianBundle, counter::Ref{Int}
)
  rows = [
    _relative_filtered_cohomology(tensor_product(F, W), counter) for
    W in reverse(_relative_wedges!(Z))
  ]
  inequalities = AffineExpr[]
  entries = long_exact_sequence_cokernel(rows, counter; inequalities)
  n = length(entries)
  system = vcat(entries, inequalities)
  for i in (dimension(Z) + 2):n
    is_zero_expr(system[i]) || _apply_equation!(system, system[i])
  end
  _propagate_intervals!(system)
  system[1:(dimension(Z) + 1)]
end

function _relative_bundle_cohomology(F::RelativeZeroLocusBundle, counter::Ref{Int})
  d = dimension(variety(F))
  isempty(F.terms) && return [AffineExpr(0) for _ in 0:d]
  row(i) = if haskey(F.terms, i)
    _relative_restriction_cohomology(variety(F), F.terms[i], counter)
  else
    [AffineExpr(0) for _ in 0:d]
  end
  entries, inequalities = _cohomology_from_presentation(row, keys(F.terms), counter)
  system = vcat(entries, inequalities)
  equation = _alternating_sum(entries) - AffineExpr(euler_characteristic(F))
  is_zero_expr(equation) || _apply_equation!(system, equation)
  _propagate_intervals!(system)
  system[1:(d + 1)]
end

"""Cohomology via Koszul and ambient presentations; unknown map ranks stay symbolic."""
function cohomology(F::RelativeZeroLocusBundle)
  entries = _relative_bundle_cohomology(F, Ref(0))
  _renumber_variables!(entries)
  Cohomology{AffineExpr}(entries, dimension(variety(F)))
end
cohomology(Z::RelativeZeroLocus, F::RelativeGrassmannianBundle) = cohomology(restrict(Z, F))

"""Hodge numbers of a smooth relative zero locus, possibly symbolic."""
function hodge_numbers(Z::RelativeZeroLocus)
  d = dimension(Z)
  M = [AffineExpr(0) for _ in 0:d, _ in 0:d]
  counter = Ref(0)
  chis = BigInt[]
  Ω = cotangent_bundle(Z)
  for p in 0:(d ÷ 2)
    F = exterior_power(Ω, p)
    M[p + 1, :] = _relative_bundle_cohomology(F, counter)
    push!(chis, euler_characteristic(F))
  end
  _propagate_hodge_constraints!(M, chis, d)
  for p in (d ÷ 2 + 1):d, q in 0:d
    M[p + 1, q + 1] = M[d - p + 1, d - q + 1]
  end
  _renumber_variables!(M)
  M
end

"""
    hochschild_cohomology(Z::RelativeZeroLocus)

Polyvector cohomology through the relative Koszul and tangent sequences.
The Hochschild--Kostant--Rosenberg row sums retain undetermined differential
ranks; symbolic variables are shared across the whole parallelogram.
"""
function hochschild_cohomology(Z::RelativeZeroLocus)
  d = dimension(Z)
  T = tangent_bundle(Z)
  K = canonical_bundle(Z)
  counter = Ref(0)
  data = Matrix{AffineExpr}(undef, d + 1, d + 1)
  serre = similar(data)
  for p in 0:d
    F = exterior_power(T, p)
    data[p + 1, :] = _relative_bundle_cohomology(F, counter)
    serre[p + 1, :] = reverse(
      _relative_bundle_cohomology(tensor_product(dual(F), K), counter)
    )
  end
  system = cat(data, serre; dims=3)
  changed = true
  while changed
    changed = false
    for p in 1:(d + 1), q in 1:(d + 1)
      equation = system[p, q, 1] - system[p, q, 2]
      is_zero_expr(equation) || (changed |= _apply_equation!(system, equation))
    end
    _propagate_intervals!(system)
  end
  data = system[:, :, 1]
  _renumber_variables!(data)
  PolyvectorParallelogram(data, d)
end

"""Topological Euler characteristic of a smooth relative zero locus."""
function euler_characteristic(Z::RelativeZeroLocus)
  Ω = cotangent_bundle(Z)
  sum(
    ((-1)^p * euler_characteristic(exterior_power(Ω, p)) for p in 0:dimension(Z));
    init=BigInt(0),
  )
end

"""Top self-intersection of a line bundle, by finite differences of its Euler characteristic."""
function degree(L::Union{ZeroLocusBundle,RelativeZeroLocusBundle})
  rank(L) == 1 || throw(ArgumentError("Degree requires a line bundle."))
  d = dimension(variety(L))
  sum(
    (
      (-1)^(d - i) * binomial(BigInt(d), i) * euler_characteristic(symmetric_power(L, i))
      for i in 0:d
    );
    init=BigInt(0),
  )
end
