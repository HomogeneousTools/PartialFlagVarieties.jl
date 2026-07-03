# ═══════════════════════════════════════════════════════════════════════════════
#  Cohomology — sheaf cohomology of equivariant bundles on G/P
#
#  Uses the Borel–Weil–Bott theorem to compute sheaf cohomology.
#  The result is a parametric Cohomology{T} where:
#    T = WeylCharacter  →  character-valued (full representation info)
#    T = BigInt          →  dimension-valued (numerical)
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
const _BWB_WEIGHT_CACHE = let b = _default_cache_budget()
  LRU{WeightLatticeElem,Union{Nothing,Tuple{Int,WeightLatticeElem}}}(;
    maxsize=_cache_maxsize(b, _DEFAULT_BWB_FRAC * 0.2),
    by=Base.summarysize,
  )
end

function _borel_weil_bott_generic(@nospecialize(λ::WeightLatticeElem))
  get!(_BWB_WEIGHT_CACHE, λ) do
    DT = typeof(λ).parameters[1]
    ρ = weyl_vector(DT)
    μ = λ + ρ
    μ_dom, d = conjugate_dominant_weight_with_length(μ)
    any(==(0), μ_dom.vec) ? nothing : (d, μ_dom - ρ)
  end
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
exactly when one Dynkin coordinate vanishes. On ``B_2`` the weight ``-ω_2`` has
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
  _borel_weil_bott_generic(λ)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition (0-based indexing)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Cohomology{T}

Sheaf cohomology of an equivariant bundle on a partial flag variety.

Entries are 0-indexed: `H[i]` returns ``\\mathrm{H}^i(G/P, E)``.

# Type parameter
- `T = WeylCharacter{DT,R}`: entries are virtual characters of ``G``
- `T = BigInt`: entries are dimensions

The object behaves like a vector indexed by cohomological degree rather than by
Julia's usual `1:length(H)` convention.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = structure_sheaf(X);

julia> H = cohomology(E);

julia> H[0]  # H⁰(ℙ⁴, 𝒪) = k
A4(0, 0, 0, 0)

julia> H[1]  # H¹(ℙ⁴, 𝒪) = 0
0
```
"""
struct Cohomology{T}
  entries::Vector{T}   # entries[i+1] = Hⁱ, i = 0, ..., dim
  dim_variety::Int      # dimension of the variety
end

# ─── 0-based indexing ────────────────────────────────────────────────────────

"""
    getindex(H::Cohomology, i::Int) -> T

Return ``\\mathrm{H}^i(G/P, E)``. Uses 0-based indexing.
"""
function Base.getindex(H::Cohomology{T}, i::Int) where {T}
  0 <= i <= H.dim_variety || throw(BoundsError(H, i))
  return H.entries[i + 1]
end

Base.length(H::Cohomology) = H.dim_variety + 1
Base.firstindex(::Cohomology) = 0
Base.lastindex(H::Cohomology) = H.dim_variety

function Base.iterate(H::Cohomology, state=0)
  state > H.dim_variety && return nothing
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
    cohomology(E::CompletelyReducibleBundle) -> Cohomology{WeylCharacter}

Compute the sheaf cohomology ``\\mathrm{H}^\\bullet(G/P, E)`` using the Borel–Weil–Bott theorem.

Returns character-valued cohomology: each ``\\mathrm{H}^i`` is a virtual character
(Weyl character) of the ambient group ``G``.

The partial flag variety is inferred from `E`.
The result is a 0-indexed [`Cohomology`](@ref) object.

# Algorithm
For each irreducible Levi component ``V_\\lambda`` of ``E``:
1. Convert to ambient weight ``\\lambda``
2. Apply BWB: `borel_weil_bott(λ)` → `(d, μ)` or `nothing`
3. If non-singular, add ``V_\\mu`` to ``\\mathrm{H}^d``

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> H = cohomology(line_bundle(X, 1));

julia> degree(H[0])  # H⁰(ℙ⁴, 𝒪(1)) = V(ω₁) of dim 5
5
```
"""
function cohomology(E::CompletelyReducibleBundle)
  d = dimension(E.variety)
  comps = components(E)
  length(comps) == 1 && return _cohomology_single(d, p_dominant_weight(only(comps)))
  sentinel = if isempty(comps)
    WeightLatticeElem(dynkin_type(marked_dynkin_type(variety(E))))
  else
    p_dominant_weight(first(comps))
  end
  _cohomology_generic(comps, d, sentinel)
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

julia> H = cohomology(line_bundle(X, 1));

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
  Cohomology{BigInt}(entries, H.dim_variety)
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

"""
    dimensions(E::CompletelyReducibleBundle) -> Cohomology{BigInt}

Compute dimension-valued cohomology directly, without first materializing the
character-valued result.

This is usually the most convenient entry point when only the dimensions of the
cohomology groups matter. As for all [`Cohomology`](@ref) objects, the result is
0-indexed: `H[0]` means ``\\mathrm{H}^0``.
"""
function dimensions(E::CompletelyReducibleBundle)
  d = dimension(E.variety)
  comps = components(E)
  isempty(comps) && return Cohomology{BigInt}(zeros(BigInt, d + 1), d)
  length(comps) == 1 && return _dimensions_single(d, p_dominant_weight(only(comps)))
  _dimensions_generic(comps, d, p_dominant_weight(first(comps)))
end

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

Compute the Euler characteristic ``\\chi(E) = \\sum_i (-1)^i \\dim \\mathrm{H}^i(G/P, E)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> H = dimensions(line_bundle(X, 1));

julia> euler_characteristic(H)
5
```
"""
function euler_characteristic(H::Cohomology{BigInt})
  sum(((-1)^i * H[i] for i in 0:(H.dim_variety)); init=BigInt(0))
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

Compute the Euler characteristic ``\\chi(G/P, E)`` directly via BWB,
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

"""
    hilbert_polynomial(E::CompletelyReducibleBundle,
                       L::CompletelyReducibleBundle;
                       max_degree::Int=20) -> Vector{Rational{BigInt}}

Compute the Hilbert polynomial of a bundle ``E`` with respect to the
polarization ``L``:
```math
P(t) = \\chi(E \\otimes L^{\\otimes t}).
```

The polarization ``L`` must be a line bundle on the same variety as ``E``;
for the result to be a genuine Hilbert polynomial, ``L`` should be ample.

Returns polynomial coefficients ``[a_0, a_1, \\ldots, a_d]`` such that
``P(t) = a_0 + a_1 t + \\cdots + a_d t^d``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(1) * projective_space(1);

julia> hilbert_polynomial(structure_sheaf(X), line_bundle(X, [1, 1]))
3-element Vector{Rational{BigInt}}:
 1
 2
 1
```
"""
function hilbert_polynomial(
  E::CompletelyReducibleBundle, L::CompletelyReducibleBundle; max_degree::Int=20
)
  variety(E) == variety(L) || throw(ArgumentError(
    "Polarization L must live on the same variety as E."
  ))
  rank_bundle(L) == 1 || throw(
    ArgumentError(
      "Polarization L must be a line bundle (rank 1), got rank $(rank_bundle(L))."
    ),
  )

  d = dimension(variety(E))
  n_points = min(d + rank_bundle(E) + 5, max_degree + 1)
  values = Rational{BigInt}[Rational{BigInt}(euler_characteristic(E))]
  Lt = L
  for _ in 1:(n_points - 1)
    push!(values, Rational{BigInt}(euler_characteristic(tensor_product(E, Lt))))
    Lt = tensor_product(Lt, L)
  end
  _newton_interpolation(values)
end

"""
    hilbert_polynomial(E::CompletelyReducibleBundle;
                       max_degree::Int=20) -> Vector{Rational{BigInt}}

Compute the Hilbert polynomial of a bundle ``\\mathcal{E}`` with respect to
the minimal ample polarization ``\\mathcal{O}(1, \\ldots, 1)`` (degree one at
every marked node); for Picard rank 1 this is the ample generator of
``\\operatorname{Pic}(X)``.  Pass an explicit line bundle as a second
argument for any other polarization.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> hilbert_polynomial(structure_sheaf(X))
4-element Vector{Rational{BigInt}}:
   1
 11//6
   1
  1//6
```
"""
function hilbert_polynomial(E::CompletelyReducibleBundle; max_degree::Int=20)
  X = variety(E)
  polarization = line_bundle(X, ones(Int, length(marked_nodes(X))))
  hilbert_polynomial(E, polarization; max_degree=max_degree)
end

"""
    _newton_interpolation(values::Vector{Rational{BigInt}}) -> Vector{Rational{BigInt}}

Given `values[i] = P(i-1)` for `i = 1, ..., n`, recover the polynomial
coefficients `[a_0, a_1, ..., a_{n-1}]` with `P(t) = Σ aₖ tᵏ` via Newton's
forward differences: `P(t) = Σ_k C(t, k) Δᵏf(0)`.
"""
function _newton_interpolation(values::Vector{Rational{BigInt}})
  n = length(values)

  deltas = Rational{BigInt}[]
  forward_differences = values
  while !isempty(forward_differences)
    push!(deltas, forward_differences[1])
    forward_differences = diff(forward_differences)
  end

  # Convert from the Newton basis C(t, k) to the monomial basis tᵏ.
  result = zeros(Rational{BigInt}, n)
  for k in 0:(n - 1)
    for (d, c) in enumerate(_binomial_poly_coeffs(k))
      result[d] += deltas[k + 1] * c
    end
  end

  # Trim trailing zeros.
  while length(result) > 1 && result[end] == 0
    pop!(result)
  end
  result
end

"""
    _binomial_poly_coeffs(k::Int) -> Vector{Rational{BigInt}}

Coefficients of C(t, k) = t(t-1)...(t-k+1)/k! as polynomial in t.
Returns [a_0, a_1, ..., a_k].
"""
function _binomial_poly_coeffs(k::Int)
  k == 0 && return Rational{BigInt}[1]

  # Start with polynomial 1
  poly = Rational{BigInt}[1]

  for j in 0:(k - 1)
    # Multiply by (t - j)
    new_poly = zeros(Rational{BigInt}, length(poly) + 1)
    for (d, c) in enumerate(poly)
      new_poly[d + 1] += c         # c * t
      new_poly[d] -= j * c         # c * (-j)
    end
    poly = new_poly
  end

  # Divide by k!
  fact = Rational{BigInt}(factorial(big(k)))
  poly ./= fact

  return poly
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Display
# ═══════════════════════════════════════════════════════════════════════════════

function Base.show(io::IO, H::Cohomology{BigInt})
  parts = String[]
  for i in 0:(H.dim_variety)
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
  for i in 0:(H.dim_variety)
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
