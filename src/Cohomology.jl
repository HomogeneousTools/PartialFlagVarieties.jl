# ═══════════════════════════════════════════════════════════════════════════════
#  Cohomology — sheaf cohomology of equivariant bundles on G/P
#
#  Uses the Borel–Weil–Bott theorem to compute sheaf cohomology.
#  The result is a parametric Cohomology{T} where:
#    T = WeylCharacter  →  character-valued (full representation info)
#    T = BigInt          →  dimension-valued (determined)
#    T = AffineExpr       →  dimension-valued (possibly symbolic)
#
#  Cohomology is stored with 0-based indexing: entry i corresponds to Hⁱ.
# ═══════════════════════════════════════════════════════════════════════════════

export Cohomology
export cohomology, dimensions
export euler_characteristic, chi
export hilbert_polynomial
export borel_weil_bott

# Names from Lie are available via the parent module.

# ═══════════════════════════════════════════════════════════════════════════════
#  Borel–Weil–Bott theorem
# ═══════════════════════════════════════════════════════════════════════════════

# Memoize the BWB kernel: every cohomology, dimension, and Euler
# characteristic call funnels through it, and the same weights recur across
# bundles (Koszul twists, plethysm summands), so the hit rate is high.
#
# The absolute statement only, for now: the key is the weight alone, which is
# sound only while every entry was folded by the same Weyl group.  Caching the
# relative fold of `borel_weil_bott(λ, nodes)` would need the node set in the
# key, since one weight has a different answer for each subsystem.
const _BWB_WEIGHT_CACHE = let b = _default_cache_budget()
  LRU{WeightLatticeElem,Union{Nothing,Tuple{Int,WeightLatticeElem}}}(;
    maxsize=_cache_maxsize(b, _DEFAULT_BWB_FRAC * 0.2),
    by=Base.summarysize,
  )
end

"""
    borel_weil_bott(λ::WeightLatticeElem{DT,R}) -> Union{Nothing, Tuple{Int, WeightLatticeElem{DT,R}}}

Apply the Borel–Weil–Bott theorem to the weight ``λ``.

Compute ``μ = λ + ρ`` and find the unique Weyl group element ``w`` such that
``w(μ)`` is dominant.  If ``μ`` is singular (any coordinate zero after
reflecting to the dominant chamber), all cohomology vanishes and
`nothing` is returned.  Otherwise return ``(d,\\, w(μ) - ρ)`` where
``d = \\ell(w)`` is the cohomological degree.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> borel_weil_bott(fundamental_weight(TypeA{2}, 1))
(0, ω1)

julia> borel_weil_bott(fundamental_weight(TypeA{2}, 1) * 0)
(0, 0)

julia> borel_weil_bott(-weyl_vector(TypeA{2})) === nothing
true
```

A *nonzero* singular weight, complementing the ``μ = 0`` example above. Note the
test is applied to the dominant representative: a *dominant* weight is singular
exactly when one Dynkin coordinate vanishes. On ``\\mathrm{B}_2`` the weight ``-ω_2`` has
``ρ``-shift ``μ = ω_1`` — dominant with one coordinate zero — so its cohomology
vanishes, while a regular weight lands in a single positive degree:

```jldoctest
julia> using PartialFlagVarieties

julia> borel_weil_bott(-fundamental_weight(TypeB{2}, 2)) === nothing  # μ = ω₁, singular
true

julia> borel_weil_bott(WeightLatticeElem(TypeB{2}, [1, -4]))          # regular: degree 2
(2, 0)
```
"""
function borel_weil_bott(λ::WeightLatticeElem{DT,R}) where {DT,R}
  # The absolute statement is the relative one with every node available, so the
  # fold lives in the two-argument method and this adds only the cache.
  get!(_BWB_WEIGHT_CACHE, λ) do
    borel_weil_bott(λ, 1:R)
  end
end

"""
    borel_weil_bott(λ::WeightLatticeElem{DT,R}, nodes) -> Union{Nothing, Tuple{Int, WeightLatticeElem{DT,R}}}

Apply the Borel–Weil–Bott theorem relative to the simple roots `nodes`, seeking
the Weyl group element in ``\\mathrm{W}_S`` rather than in ``\\mathrm{W}_{\\mathrm{G}}``.
This is the fibrewise statement used by [`pushforward`](@ref); passing every node
recovers the absolute one.

Not memoized, unlike the absolute case. Two reasons: it is called once per
irreducible summand rather than from an inner loop, and the cache behind
[`borel_weil_bott(λ)`](@ref) keys on the weight alone, so it cannot hold relative
answers as it stands. A weight has one answer per subsystem ``S``, so the node
set would have to enter the key.

Measured, memoizing this would in any case be a pessimization for its current
caller: hashing a summand costs more than the fold it saves, because
``\\mathrm{W}_S`` is small.

The shift is still ``\\rho = \\rho_{\\mathrm{G}}``: the difference
``\\rho_{\\mathrm{G}} - \\rho_S`` pairs to zero with every coroot in ``S``, so it
is ``\\mathrm{W}_S``-invariant and drops out.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> λ = WeightLatticeElem(TypeA{2}, [-2, 1]);

julia> borel_weil_bott(λ)          # regular for G: one reflection
(1, 0)

julia> borel_weil_bott(λ, (2,))    # already dominant for the second node alone
(0, -2ω1 + ω2)
```
"""
function borel_weil_bott(λ::WeightLatticeElem, nodes)
  # Deliberately not `@nospecialize`, unlike the bundle loops further down: this
  # is four lines and the compile cost sits in Semisimple's `@inline` fold, which
  # specialises regardless.  Measured, annotating it saves no latency and costs
  # throughput on both callers, the cached absolute one included.
  ρ = weyl_vector(typeof(λ).parameters[1])
  μ_dom, d = conjugate_dominant_weight_with_length(λ + ρ, nodes)
  any(s -> iszero(coefficients(μ_dom)[s]), nodes) ? nothing : (d, μ_dom - ρ)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition (0-based indexing)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Cohomology{T}

A degree-indexed container for sheaf cohomology or higher direct images.

Entries are 0-indexed: `H[i]` returns ``\\mathrm{H}^i(\\mathrm{G}/\\mathrm{P}, \\mathcal{E})``.

# Type parameter
- `T = WeylCharacter{DT,R}`: entries are virtual characters of ``\\mathrm{G}``
- `T = BigInt`: entries are dimensions
- `T = AffineExpr`: entries are dimensions that may depend on undetermined
  induced-map ranks, as for [`FilteredBundle`](@ref) and
  [`ZeroLocusBundle`](@ref) cohomology
- `T = CompletelyReducibleBundle`: entries are the higher direct images
  ``\\mathrm{R}^iq_*`` of [`pushforward`](@ref), which live on the target of a
  projection rather than being cohomology of anything. For that parametrization
  `H[i]` is a bundle and the index runs to the relative dimension.

The object behaves like a vector indexed by cohomological degree rather than by
Julia's usual `1:length(H)` convention.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = structure_sheaf(X);

julia> H = cohomology(E; characters=true);

julia> H[0]  # H⁰(ℙ⁴, 𝒪) = k
A4(0, 0, 0, 0)

julia> H[1]  # H¹(ℙ⁴, 𝒪) = 0
0
```
"""
struct Cohomology{T}
  entries::Vector{T}   # entries[i+1] = the degree-i part, i = 0, ..., max_degree
  # Highest degree that can carry something: `dim X` for sheaf cohomology on
  # `X`, the relative dimension for the higher direct images of `pushforward`.
  max_degree::Int
end

# ─── 0-based indexing ────────────────────────────────────────────────────────

"""
    getindex(H::Cohomology, i::Int) -> T

Return ``\\mathrm{H}^i(\\mathrm{G}/\\mathrm{P}, \\mathcal{E})``. Uses 0-based indexing.
"""
function Base.getindex(H::Cohomology{T}, i::Int) where {T}
  0 <= i <= H.max_degree || throw(BoundsError(H, i))
  return H.entries[i + 1]
end

Base.length(H::Cohomology) = H.max_degree + 1
Base.firstindex(::Cohomology) = 0
Base.lastindex(H::Cohomology) = H.max_degree

function Base.iterate(H::Cohomology, state=0)
  state > H.max_degree && return nothing
  return (H[state], state + 1)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Character-valued cohomology computation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _weight_counts(comps, sentinel::WeightLatticeElem) -> Dict

Count equal ambient weights in `comps` without specializing the caller on the
rank-specific `WeightLatticeElem{DT,R}` type.

The first-pass latency for workloads that iterate over many Grassmannian ranks
is dominated by recompiling the cohomology loops for each new `{DT,R}` pair.
Building the dictionary with the runtime key type `typeof(sentinel)` preserves
fast concrete hashing while keeping the outer methods generic.
"""
function _weight_counts(comps::Vector{IrrepLevi}, sentinel::WeightLatticeElem)
  Base.@nospecialize sentinel

  weight_type = typeof(sentinel)
  weight_counts = Dict{weight_type,Int}()
  for comp in comps
    λ = p_dominant_weight(comp)
    weight_counts[λ] = get(weight_counts, λ, 0) + 1
  end
  weight_counts
end

"""
    _cohomology_generic(comps, d, sentinel::WeightLatticeElem) -> Cohomology

Latency-optimized inner loop for character-valued cohomology.

The sentinel determines the runtime weight type for the per-call dictionary and
the zero character used to initialise the cohomology groups, but the method
itself is kept generic to avoid recompiling it once per ambient rank.
"""
function _cohomology_generic(
  comps::Vector{IrrepLevi}, d::Int, sentinel::WeightLatticeElem
)
  Base.@nospecialize sentinel

  DT = typeof(sentinel).parameters[1]
  entries = [WeylCharacter(DT) for _ in 0:d]
  weight_counts = _weight_counts(comps, sentinel)

  for (λ, mult) in weight_counts
    result = borel_weil_bott(λ)
    result === nothing && continue
    (deg, μ) = result
    0 <= deg <= d || continue
    χμ = WeylCharacter(μ)
    for _ in 1:mult
      add!(entries[deg + 1], χμ)
    end
  end

  Cohomology{eltype(entries)}(entries, d)
end

function _cohomology_single(d::Int, λ::WeightLatticeElem)
  Base.@nospecialize λ

  DT = typeof(λ).parameters[1]
  entries = [WeylCharacter(DT) for _ in 0:d]

  result = borel_weil_bott(λ)
  if result !== nothing
    deg, μ = result
    0 <= deg <= d && add!(entries[deg + 1], WeylCharacter(μ))
  end

  Cohomology{eltype(entries)}(entries, d)
end

"""
    cohomology(E::CompletelyReducibleBundle; characters=false) -> Cohomology

Compute the sheaf cohomology ``\\mathrm{H}^\\bullet(\\mathrm{G}/\\mathrm{P}, \\mathcal{E})`` using the Borel–Weil–Bott theorem.

By default the result is dimension-valued (`Cohomology{BigInt}`): each ``\\mathrm{H}^i``
is given by its dimension, computed directly without materializing the characters.
Pass `characters=true` to get character-valued cohomology (`Cohomology{WeylCharacter}`)
instead, where each ``\\mathrm{H}^i`` is a virtual character (Weyl character) of the
ambient group ``\\mathrm{G}``.

The partial flag variety is inferred from `E`.
The result is a 0-indexed [`Cohomology`](@ref) object.

# Algorithm
For each irreducible Levi component ``V_\\lambda`` of ``\\mathcal{E}``:
1. Convert to ambient weight ``\\lambda``
2. Apply BWB: `borel_weil_bott(λ)` → `(d, μ)` or `nothing`
3. If non-singular, add ``V_\\mu`` (or just its dimension) to ``\\mathrm{H}^d``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> cohomology(line_bundle(X, 1))[0]  # h⁰(ℙ⁴, 𝒪(1)) = 5
5

julia> degree(cohomology(line_bundle(X, 1); characters=true)[0])  # V(ω₁), dim 5
5
```
"""
function cohomology(E::CompletelyReducibleBundle; characters::Bool=false)
  d = dimension(E.variety)
  comps = components(E)
  if characters
    length(comps) == 1 && return _cohomology_single(d, p_dominant_weight(only(comps)))
    sentinel = if isempty(comps)
      WeightLatticeElem(dynkin_type(marked_dynkin_type(variety(E))))
    else
      p_dominant_weight(first(comps))
    end
    return _cohomology_generic(comps, d, sentinel)
  end
  isempty(comps) && return Cohomology{BigInt}(zeros(BigInt, d + 1), d)
  length(comps) == 1 && return _dimensions_single(d, p_dominant_weight(only(comps)))
  _dimensions_generic(comps, d, p_dominant_weight(first(comps)))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Dimension-valued cohomology
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dimensions(H::Cohomology{<:WeylCharacter}) -> Cohomology{BigInt}

Convert character-valued cohomology to dimension-valued.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> H = cohomology(line_bundle(X, 1); characters=true);

julia> dims = dimensions(H);

julia> dims[0]
5
```
"""
function dimensions(H::Cohomology{<:WeylCharacter})
  entries = BigInt[
    sum((BigInt(mult) * degree(weight) for (weight, mult) in χ); init=BigInt(0)) for
    χ in H.entries
  ]
  Cohomology{BigInt}(entries, H.max_degree)
end

"""
    _dimensions_generic(comps, d, sentinel::WeightLatticeElem) -> Cohomology{BigInt}

Latency-optimized inner loop for dimension-valued cohomology.

This keeps the method generic across ambient ranks while still using a
concrete dictionary key type determined at runtime from `sentinel`.
"""
function _dimensions_generic(
  comps::Vector{IrrepLevi}, d::Int, sentinel::WeightLatticeElem
)
  Base.@nospecialize sentinel

  entries = zeros(BigInt, d + 1)
  weight_counts = _weight_counts(comps, sentinel)
  for (λ, mult) in weight_counts
    result = borel_weil_bott(λ)
    result === nothing && continue
    deg, μ = result
    if 0 <= deg <= d
      entries[deg + 1] += mult * degree(μ)
    end
  end
  Cohomology{BigInt}(entries, d)
end

function _dimensions_single(d::Int, λ::WeightLatticeElem)
  Base.@nospecialize λ

  entries = zeros(BigInt, d + 1)
  result = borel_weil_bott(λ)
  if result !== nothing
    deg, μ = result
    if 0 <= deg <= d
      entries[deg + 1] = degree(μ)
    end
  end
  Cohomology{BigInt}(entries, d)
end

# `dimensions(E)` on a bundle is superseded by `cohomology(E)`, which is now
# dimension-valued by default. Kept as a deprecated alias for one minor version.
@deprecate dimensions(E::CompletelyReducibleBundle) cohomology(E)

# ═══════════════════════════════════════════════════════════════════════════════
#  iszero
# ═══════════════════════════════════════════════════════════════════════════════

"""
    iszero(H::Cohomology) -> Bool

Return `true` when all cohomology groups vanish.
"""
Base.iszero(H::Cohomology{BigInt}) = all(iszero, H.entries)

Base.iszero(H::Cohomology{<:WeylCharacter}) = all(χ -> isempty(χ.terms), H.entries)

# ═══════════════════════════════════════════════════════════════════════════════
#  Euler characteristic
# ═══════════════════════════════════════════════════════════════════════════════

"""
    euler_characteristic(H::Cohomology{BigInt}) -> BigInt

Compute the Euler characteristic ``\\chi(\\mathcal{E}) = \\sum_i (-1)^i \\dim \\mathrm{H}^i(\\mathrm{G}/\\mathrm{P}, \\mathcal{E})``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> H = cohomology(line_bundle(X, 1));

julia> euler_characteristic(H)
5
```
"""
function euler_characteristic(H::Cohomology{BigInt})
  sum(((-1)^i * H[i] for i in 0:(H.max_degree)); init=BigInt(0))
end

function _alternating_euler_characteristic(cohos::AbstractVector{<:Cohomology{BigInt}})
  sum(
    ((-1)^(i - 1) * euler_characteristic(H) for (i, H) in enumerate(cohos)); init=BigInt(0)
  )
end

function euler_characteristic(H::Cohomology{<:WeylCharacter})
  return euler_characteristic(dimensions(H))
end

"""
    chi(H::Cohomology) -> BigInt

Synonym for [`euler_characteristic(H)`](@ref).

!!! warning "Two different invariants"
    For `Cohomology{BigInt}` and `Cohomology{<:WeylCharacter}` this is the
    alternating sum of the dimensions of the cohomology groups. For
    `Cohomology{CompletelyReducibleBundle}` the entries are bundles, not
    cohomology, so it is the alternating sum of their *ranks*: the rank of the
    class of ``\\mathrm{R}q_*\\mathcal{E}`` in the Grothendieck group, which is
    not ``\\chi`` of any sheaf. Generic code taking an `H::Cohomology` and
    reporting `chi(H)` gets whichever of the two matches the element type.
"""
chi(H::Cohomology) = euler_characteristic(H)

"""
    _euler_characteristic_generic(comps, sentinel::WeightLatticeElem) -> BigInt

Latency-optimized Euler-characteristic loop that avoids recompiling once per
ambient rank while still deduplicating equal weights before the BWB calls.
"""
function _euler_characteristic_generic(
  comps::Vector{IrrepLevi}, sentinel::WeightLatticeElem
)::BigInt
  Base.@nospecialize sentinel

  weight_counts = _weight_counts(comps, sentinel)
  result = BigInt(0)
  for (λ, mult) in weight_counts
    bwb = borel_weil_bott(λ)
    bwb === nothing && continue
    (deg, μ) = bwb
    result += (iseven(deg) ? 1 : -1) * mult * degree(μ)
  end
  result
end

function _euler_characteristic_single(λ::WeightLatticeElem)::BigInt
  Base.@nospecialize λ

  bwb = borel_weil_bott(λ)
  bwb === nothing && return BigInt(0)
  deg, μ = bwb
  (iseven(deg) ? 1 : -1) * degree(μ)
end

"""
    euler_characteristic(E::CompletelyReducibleBundle) -> BigInt

Compute the Euler characteristic ``\\chi(\\mathrm{G}/\\mathrm{P}, \\mathcal{E})`` directly via BWB,
without constructing intermediate `WeylCharacter` objects.

For each irreducible Levi component, converts to the ambient weight,
applies the Borel–Weil–Bott theorem, and accumulates
``(-1)^{\\ell(w)} \\dim V_{w(\\lambda+\\rho)-\\rho}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> euler_characteristic(structure_sheaf(X))
1
```
"""
function euler_characteristic(E::CompletelyReducibleBundle)
  comps = components(E)
  isempty(comps) && return BigInt(0)
  length(comps) == 1 && return _euler_characteristic_single(p_dominant_weight(only(comps)))
  _euler_characteristic_generic(comps, p_dominant_weight(first(comps)))
end

"""
    chi(E::CompletelyReducibleBundle) -> BigInt

Synonym for [`euler_characteristic(E)`](@ref).
"""
chi(E::CompletelyReducibleBundle) = euler_characteristic(E)

# ═══════════════════════════════════════════════════════════════════════════════
#  Hilbert polynomial
# ═══════════════════════════════════════════════════════════════════════════════

# Require a polarization to be a line bundle on the given ambient variety.
function _check_polarization(X::PartialFlagVariety, L::CompletelyReducibleBundle)
  variety(L) == X ||
    throw(ArgumentError("Polarization L must live on the same ambient variety."))
  rank(L) == 1 || throw(ArgumentError("Polarization L must be a line bundle (rank 1)."))
  nothing
end

"""
    hilbert_polynomial(E::CompletelyReducibleBundle[, L::CompletelyReducibleBundle])
      -> Polynomial

Compute the Hilbert polynomial of a bundle ``\\mathcal{E}`` with respect to the
polarization ``\\mathcal{L}``:
```math
\\begin{equation}
P(t) = \\chi(\\mathcal{E} \\otimes \\mathcal{L}^{\\otimes t}).
\\end{equation}
```

The polarization ``\\mathcal{L}`` must be a line bundle on the same variety as ``\\mathcal{E}``;
for the result to be a genuine Hilbert polynomial, ``\\mathcal{L}`` should be ample.
Omitting `L` uses the minimal ample polarization ``\\mathcal{O}(1, \\ldots, 1)``,
with degree one at every marked node.
For Picard rank one, this is the ample generator of ``\\operatorname{Pic}(X)``.

For an irreducible Levi summand with ambient weight ``\\lambda`` and polarization
weight ``\\mu``, Borel–Weil–Bott and the Weyl dimension formula give
```math
\\begin{equation}
\\chi(\\mathcal{E}_\\lambda \\otimes \\mathcal{L}^{\\otimes t})
= \\prod_{\\alpha^\\vee > 0}
\\frac{\\langle \\lambda + \\rho + t\\mu, \\alpha^\\vee \\rangle}
     {\\langle \\rho, \\alpha^\\vee \\rangle}.
\\end{equation}
```
Here ``\\rho`` is the Weyl vector and the product runs over all positive coroots.
The product vanishes at singular weights; at regular weights its Weyl sign
is exactly the cohomological sign in Borel–Weil–Bott.
Thus it computes the Euler characteristic for every integer twist, including
negative twists, with exact rational coefficients.
For ample `L`, Serre vanishing identifies it with the usual Hilbert polynomial.
See Fulton–Harris, *Representation Theory: A First Course*,
[§§23.3 and 24.1](https://doi.org/10.1007/978-1-4612-0979-9),
and the [mathematical background](@ref hilbert-polynomials).

Returns a [`Polynomial`](@ref), evaluated with `P(t)`.
Use [`constant_coefficient`](@ref), [`leading_coefficient`](@ref), and
[`coefficient`](@ref) for individual coefficients, or
[`coefficients(::Polynomial)`](@ref) for the coefficient vector in increasing degree order.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(1) * projective_space(1);

julia> hilbert_polynomial(structure_sheaf(X), line_bundle(X, [1, 1]))
t^2 + 2*t + 1

julia> hilbert_polynomial(structure_sheaf(projective_space(3)))
1//6*t^3 + t^2 + 11//6*t + 1
```
"""
function hilbert_polynomial(
  E::CompletelyReducibleBundle, L::CompletelyReducibleBundle
)
  _check_polarization(variety(E), L)
  result = zero(Polynomial)
  coroots = positive_coroots(RootSystem(dynkin_type(variety(E))))
  μ = p_dominant_weight(only(components(L))) # Weight of the polarization L.
  for (rep, mult) in _to_counts(E)
    λ = p_dominant_weight(rep)
    polynomial = one(Polynomial)
    for coroot in coroots
      # Weights use fundamental-weight coordinates and coroots use
      # simple-coroot coordinates, so their natural pairing is the coordinate
      # dot product. Semisimple.dot instead computes a metric inner product.
      # Promote before arithmetic so the Weyl factors cannot overflow Int.
      coroot_coefficients = coefficients(coroot)
      indices = eachindex(coroot_coefficients)

      # Since ρ has all fundamental-weight coordinates equal to one,
      # these are <ρ, α∨>, <λ + ρ, α∨>, and <μ, α∨>, respectively.
      denominator = sum(BigInt, coroot_coefficients)
      constant = sum((BigInt(λ[i]) + 1) * coroot_coefficients[i] for i in indices)
      slope = sum(BigInt(μ[i]) * coroot_coefficients[i] for i in indices)

      # Multiply by (<λ + ρ, α∨> + t<μ, α∨>) / <ρ, α∨>.
      # A zero factor means every twist is singular and hence acyclic.
      polynomial *= Polynomial([constant, slope]) / denominator
      iszero(polynomial) && break
    end
    result += mult * polynomial
  end
  result
end

"""
    hilbert_polynomial(F::FilteredBundle[, L::CompletelyReducibleBundle])
      -> Polynomial

Compute the Hilbert polynomial by adding the symbolic Weyl polynomials of
its associated graded pieces. Additivity makes the result independent of
extension maps. The default polarization has degree one at every marked node.
"""
function hilbert_polynomial(F::FilteredBundle, L::CompletelyReducibleBundle)
  _check_polarization(variety(F), L)
  sum(piece -> hilbert_polynomial(piece, L), graded_pieces(F); init=zero(Polynomial))
end

hilbert_polynomial(F::FilteredBundle) =
  hilbert_polynomial(F, line_bundle(variety(F), ones(Int, picard_rank(variety(F)))))

function hilbert_polynomial(E::CompletelyReducibleBundle)
  X = variety(E)
  L = line_bundle(X, ones(Int, picard_rank(X)))
  hilbert_polynomial(E, L)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Display
# ═══════════════════════════════════════════════════════════════════════════════

function Base.show(io::IO, H::Cohomology{BigInt})
  parts = String[]
  for i in 0:(H.max_degree)
    v = H[i]
    v == 0 && continue
    push!(parts, "H$(_superscript(i)) = $v")
  end
  if isempty(parts)
    print(io, "H* = 0")
  else
    print(io, join(parts, "\n"))
  end
end

function Base.show(io::IO, H::Cohomology{<:WeylCharacter})
  parts = String[]
  for i in 0:(H.max_degree)
    v = H[i]
    isempty(v.terms) && continue
    push!(parts, "H$(_superscript(i)) = $(sprint(show, v))")
  end
  if isempty(parts)
    print(io, "H* = 0")
  else
    print(io, join(parts, "\n"))
  end
end

const _SUPERSCRIPT_DIGITS = Dict(
  '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
  '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
)

_superscript(n::Int) = String([get(_SUPERSCRIPT_DIGITS, c, c) for c in string(n)])
