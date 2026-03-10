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

"""
    borel_weil_bott(λ::WeightLatticeElem{DT,R}) -> Union{Nothing, Tuple{Int, WeightLatticeElem{DT,R}}}

Apply the Borel–Weil–Bott theorem to the weight `λ`.

Compute `μ = λ + ρ` and find the unique Weyl group element `w` such that
`w(μ)` is dominant.  If `μ` is singular (any coordinate zero after
reflecting to the dominant chamber), all cohomology vanishes and
`nothing` is returned.  Otherwise return `(d, w(μ) - ρ)` where
`d = ℓ(w)` is the cohomological degree.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> borel_weil_bott(fundamental_weight(TypeA{2}, 1))
(0, ω1)

julia> borel_weil_bott(fundamental_weight(TypeA{2}, 1) * 0)
(0, 0)

julia> borel_weil_bott(-Lie.weyl_vector(TypeA{2})) === nothing
true
```
"""
function borel_weil_bott(λ::WeightLatticeElem{DT,R}) where {DT,R}
  ρ = Lie.weyl_vector(DT)
  μ = λ + ρ
  μ_dom, d = Lie.conjugate_dominant_weight_with_length(μ)
  any(==(0), μ_dom.vec) && return nothing
  return (d, μ_dom - ρ)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition (0-based indexing)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Cohomology{T}

Sheaf cohomology of an equivariant bundle on a partial flag variety.

Entries are 0-indexed: `H[i]` returns ``H^i(G/P, E)``.

# Type parameter
- `T = WeylCharacter{DT,R}`: entries are virtual characters of ``G``
- `T = BigInt`: entries are dimensions

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

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

Return ``H^i(G/P, E)``. Uses 0-based indexing.
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
    cohomology(E::CompletelyReducibleBundle) -> Cohomology{WeylCharacter}

Compute the sheaf cohomology ``H^*(G/P, E)`` using the Borel–Weil–Bott theorem.

Returns character-valued cohomology: each ``H^i`` is a virtual character
(Weyl character) of the ambient group ``G``.

The partial flag variety is inferred from `E`.

# Algorithm
For each irreducible Levi component ``V_\\lambda`` of ``E``:
1. Convert to ambient weight ``\\lambda``
2. Apply BWB: `borel_weil_bott(λ)` → `(d, μ)` or `nothing`
3. If non-singular, add ``V_\\mu`` to ``H^d``

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> H = cohomology(line_bundle(X, 1));

julia> degree(H[0])  # H⁰(ℙ⁴, 𝒪(1)) = V(ω₁) of dim 5
5
```
"""
function cohomology(E::CompletelyReducibleBundle)
  mdt = marked_dynkin_type(variety(E))
  DT = dynkin_type(mdt)
  R = rank(DT)
  d = dimension(E.variety)

  # Initialize cohomology entries as zero characters
  entries = [WeylCharacter(DT) for _ in 0:d]

  weight_counts = Dict{WeightLatticeElem,Int}()
  for comp in components(E)
    λ = to_ambient_weight(mdt, comp)
    weight_counts[λ] = get(weight_counts, λ, 0) + 1
  end

  for (λ, mult) in weight_counts
    result = borel_weil_bott(λ)

    if result !== nothing
      (deg, μ) = result
      if 0 <= deg <= d
        χμ = WeylCharacter(μ)
        for _ in 1:mult
          add!(entries[deg + 1], χμ)
        end
      end
    end
  end

  return Cohomology{WeylCharacter{DT,R}}(entries, d)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Dimension-valued cohomology
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dimensions(H::Cohomology{<:WeylCharacter}) -> Cohomology{BigInt}

Convert character-valued cohomology to dimension-valued.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> H = cohomology(line_bundle(X, 1));

julia> dims = dimensions(H);

julia> dims[0]
5
```
"""
function dimensions(H::Cohomology{<:WeylCharacter{DT,R}}) where {DT,R}
  dim_entries = BigInt[]
  for χ in H.entries
    if isempty(χ.terms)
      push!(dim_entries, BigInt(0))
    else
      d = BigInt(0)
      for (hw, mult) in χ
        d += mult * degree(hw)
      end
      push!(dim_entries, d)
    end
  end
  return Cohomology{BigInt}(dim_entries, H.dim_variety)
end

"""
    dimensions(E::CompletelyReducibleBundle) -> Cohomology{BigInt}

Compute dimension-valued cohomology directly.
"""
function dimensions(E::CompletelyReducibleBundle)
  mdt = marked_dynkin_type(variety(E))
  d = dimension(E.variety)
  entries = zeros(BigInt, d + 1)

  weight_counts = Dict{WeightLatticeElem,Int}()
  for comp in components(E)
    λ = to_ambient_weight(mdt, comp)
    weight_counts[λ] = get(weight_counts, λ, 0) + 1
  end

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

# ═══════════════════════════════════════════════════════════════════════════════
#  iszero
# ═══════════════════════════════════════════════════════════════════════════════

"""
    iszero(H::Cohomology) -> Bool

Return `true` when all cohomology groups vanish.
"""
function Base.iszero(H::Cohomology{BigInt})
  return all(==(BigInt(0)), H.entries)
end

function Base.iszero(H::Cohomology{<:WeylCharacter})
  return all(isempty(v.terms) for v in H.entries)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Euler characteristic
# ═══════════════════════════════════════════════════════════════════════════════

"""
    euler_characteristic(H::Cohomology{BigInt}) -> BigInt

Compute the Euler characteristic ``\\chi(E) = \\sum_i (-1)^i \\dim H^i(G/P, E)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> H = dimensions(line_bundle(X, 1));

julia> euler_characteristic(H)
5
```
"""
function euler_characteristic(H::Cohomology{BigInt})
  result = BigInt(0)
  for i in 0:H.dim_variety
    result += (iseven(i) ? 1 : -1) * H[i]
  end
  return result
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
    euler_characteristic(E::CompletelyReducibleBundle) -> BigInt

Compute the Euler characteristic ``\\chi(G/P, E)`` directly via BWB,
without constructing intermediate `WeylCharacter` objects.

For each irreducible Levi component, converts to the ambient weight,
applies the Borel–Weil–Bott theorem, and accumulates
``(-1)^{\\ell(w)} \\dim V_{w(\\lambda+\\rho)-\\rho}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> euler_characteristic(structure_sheaf(X))
1
```
"""
function euler_characteristic(E::CompletelyReducibleBundle)
  mdt = marked_dynkin_type(variety(E))
  # Deduplicate: count multiplicities of identical ambient weights
  weight_counts = Dict{WeightLatticeElem,Int}()
  for comp in components(E)
    λ = to_ambient_weight(mdt, comp)
    weight_counts[λ] = get(weight_counts, λ, 0) + 1
  end
  result = BigInt(0)
  for (λ, mult) in weight_counts
    bwb = borel_weil_bott(λ)
    bwb === nothing && continue
    (deg, μ) = bwb
    result += (iseven(deg) ? 1 : -1) * mult * degree(μ)
  end
  result
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
    hilbert_polynomial(E::CompletelyReducibleBundle;
                       max_degree::Int=20) -> Vector{Rational{BigInt}}

Compute the Hilbert polynomial of a bundle ``E`` on a generalized Grassmannian.

For ``\\mathcal{O}(1)`` the Serre twist ``\\chi(E(t)) = \\sum_i (-1)^i \\dim H^i(E \\otimes \\mathcal{L}(t))``.

Returns polynomial coefficients ``[a_0, a_1, ..., a_d]`` such that
``P(t) = a_0 + a_1 t + ... + a_d t^d``.

Requires the variety to be a generalized Grassmannian (one marked node).

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));  # ℙ⁴

julia> coeffs = hilbert_polynomial(structure_sheaf(X));

julia> length(coeffs) > 0
true
```
"""
function hilbert_polynomial(E::CompletelyReducibleBundle;
  max_degree::Int=20,
)
  marked = marked_nodes(variety(E))
  length(marked) == 1 || throw(
    ArgumentError(
      "Hilbert polynomial requires a generalized Grassmannian (1 marked node)"
    ),
  )

  d = dimension(E.variety)
  # By Riemann-Roch-Hirzebruch or direct computation,
  # evaluate χ(E(t)) at enough integer points to interpolate

  n_points = min(d + rank_bundle(E) + 5, max_degree + 1)
  values = Rational{BigInt}[]

  for t in 0:(n_points - 1)
    Et = twist(E, 1, t)
    push!(values, Rational{BigInt}(euler_characteristic(Et)))
  end

  # Lagrange interpolation to recover polynomial
  coeffs = _lagrange_interpolation(values)
  return coeffs
end

"""
    _lagrange_interpolation(values::Vector{Rational{BigInt}}) -> Vector{Rational{BigInt}}

Given `values[i] = P(i-1)` for `i = 1, ..., n`, recover the polynomial coefficients
`[a_0, a_1, ..., a_{n-1}]` such that `P(t) = Σ aₖ tᵏ`.
"""
function _lagrange_interpolation(values::Vector{Rational{BigInt}})
  n = length(values)
  # Newton's forward differences
  # Δ⁰f(0) = f(0), Δ¹f(0) = f(1) - f(0), etc.
  # P(t) = Σ_{k=0}^{n-1} C(t, k) * Δᵏf(0)
  # where C(t, k) = t*(t-1)*...*(t-k+1)/k!

  diffs = copy(values)
  deltas = Rational{BigInt}[diffs[1]]

  for order in 1:(n - 1)
    new_diffs = Rational{BigInt}[]
    for i in 1:(length(diffs) - 1)
      push!(new_diffs, diffs[i + 1] - diffs[i])
    end
    push!(deltas, new_diffs[1])
    diffs = new_diffs
  end

  # Convert from Newton basis C(t, k) to monomial basis tᵏ
  # C(t, k) = (1/k!) tᵏ + lower terms
  # Use symbolic expansion
  # coeffs[d+1] = coefficient of t^d
  result = zeros(Rational{BigInt}, n)
  for k in 0:(n - 1)
    # Expand C(t, k) * deltas[k+1]
    # C(t, k) = Σ_{d=0}^{k} s(k, d) / k! * t^d  where s are Stirling
    # Simpler: compute the polynomial t(t-1)...(t-k+1)/k! directly
    binom_poly = _binomial_poly_coeffs(k)
    for (d, c) in enumerate(binom_poly)
      result[d] += deltas[k + 1] * c
    end
  end

  # Trim trailing zeros
  while length(result) > 1 && result[end] == 0
    pop!(result)
  end

  return result
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
  for i in 0:H.dim_variety
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
  for i in 0:H.dim_variety
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

function _superscript(n::Int)
  digits = Dict(
    '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
  )
  return String([get(digits, c, c) for c in string(n)])
end
