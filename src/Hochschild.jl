# ═══════════════════════════════════════════════════════════════════════════════
#  Hochschild — Hodge numbers, twisted Hodge numbers, and Hochschild cohomology
#
#  For smooth projective varieties:
#  - hodge_numbers(X): the Hodge diamond h^{p,q}(X)
#  - twisted_hodge_numbers(X, j): twisted Hodge numbers h^q(X, Ω^p(j))
#  - hochschild_cohomology(X): HKR decomposition HH^n = ⊕_{p+q=n} H^q(∧^p T)
#
#  For rational homogeneous varieties G/P:
#  - h^{p,q} = 0 for p ≠ q (rationality)
#  - h^{p,p} = b_{2p} (Betti numbers)
# ═══════════════════════════════════════════════════════════════════════════════

export hodge_numbers, twisted_hodge_numbers
export hochschild_cohomology
export PolyvectorParallelogram

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge numbers  h^{p,q}(X)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers(X::PartialFlagVariety) -> Matrix{BigInt}

Compute the Hodge numbers ``h^{p,q}(X) = \\dim H^q(X, \\Omega^p_X)``
of the partial flag variety ``X = G/P``.

Since ``G/P`` is rational (and simply connected), the Hodge diamond is
diagonal: ``h^{p,q} = 0`` for ``p \\neq q``, and ``h^{p,p} = b_{2p}``
where ``b_i`` are the Betti numbers.

Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = h^{p,q}``,
with ``d = \\dim X``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 4);

julia> H = hodge_numbers(X);

julia> H[1, 1]  # h^{0,0}
1

julia> H[2, 2]  # h^{1,1}
1

julia> H[1, 2]  # h^{0,1}
0
```
"""
function hodge_numbers(X::PartialFlagVariety)
  d = dimension(X)
  betti = betti_numbers(X)

  # For G/P: h^{p,q} = δ_{p,q} * b_{2p}
  # betti_numbers returns only even Betti numbers: betti[i] = b_{2(i-1)}
  H = zeros(BigInt, d + 1, d + 1)
  for p in 0:d
    if p + 1 <= length(betti)
      H[p + 1, p + 1] = betti[p + 1]  # betti[p+1] = b_{2p}
    end
  end
  H
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Twisted Hodge numbers  h^q(X, Ω^p(j))
# ═══════════════════════════════════════════════════════════════════════════════

"""
    twisted_hodge_numbers(X::PartialFlagVariety, j::Integer) -> Matrix{BigInt}

Compute the twisted Hodge numbers ``h^q(X, \\Omega^p_X(j))`` for the
twist ``j``.

Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = h^q(X, \\Omega^p(j))``.

Requires ``X`` to be a generalized Grassmannian (Picard rank 1).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(3);

julia> H = twisted_hodge_numbers(X, 0);

julia> H[1, 1]  # h^0(ℙ³, 𝒪)
1
```
"""
function twisted_hodge_numbers(X::PartialFlagVariety{MDT}, j::Integer) where {
  MDT<:MarkedDynkinType{DT,Marked}
} where {DT,Marked}
  length(Marked) == 1 || throw(ArgumentError(
    "twisted_hodge_numbers requires Picard rank 1"
  ))

  d = dimension(X)
  H = zeros(BigInt, d + 1, d + 1)

  for p in 0:d
    # Ω^p(j) = ∧^p Ω ⊗ O(j)
    Omega_p = exterior_power(cotangent_bundle(X), p)
    Omega_p_j = twist(Omega_p, 1, Int(j))
    coh = dimensions(Omega_p_j)
    for q in 0:d
      H[p + 1, q + 1] = coh[q]
    end
  end

  # Apply symmetries and vanishings:
  # Serre duality: h^q(Ω^p(j)) = h^{d-q}(Ω^{d-p}(-j+c₁))
  # These are already encoded in the BWB computation.

  H
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Polyvector Parallelogram — display type for Hochschild cohomology
# ═══════════════════════════════════════════════════════════════════════════════

"""
    PolyvectorParallelogram

The Hochschild–Kostant–Rosenberg decomposition of Hochschild cohomology:

``\\mathrm{HH}^n(X) = \\bigoplus_{p+q=n} H^q(X, \\bigwedge^p T_X)``

Stored as a matrix `data[p+1, q+1] = h^q(X, ∧^p T_X)`.

# Fields
- `data::Matrix{BigInt}`: the HKR decomposition matrix
- `dim::Int`: the dimension of the variety
"""
struct PolyvectorParallelogram
  data::Matrix{BigInt}  # data[p+1, q+1] = h^q(X, ∧^p T_X)
  dim::Int
end

"""
    getindex(P::PolyvectorParallelogram, p::Int, q::Int) -> BigInt

Return ``h^q(X, \\bigwedge^p T_X)``. Uses 0-based indexing.
"""
function Base.getindex(P::PolyvectorParallelogram, p::Int, q::Int)
  0 <= p <= P.dim || return BigInt(0)
  0 <= q <= P.dim || return BigInt(0)
  P.data[p + 1, q + 1]
end

"""
    euler_characteristic(P::PolyvectorParallelogram) -> BigInt

Compute the Euler characteristic of Hochschild cohomology:
``\\chi(\\mathrm{HH}^*) = \\sum_{p,q} (-1)^{p+q} h^q(X, \\bigwedge^p T_X)``
"""
function euler_characteristic(P::PolyvectorParallelogram)
  result = BigInt(0)
  for p in 0:P.dim
    for q in 0:P.dim
      result += (-1)^(p + q) * P[p, q]
    end
  end
  result
end

"""
    hochschild_dimension(P::PolyvectorParallelogram, n::Int) -> BigInt

Compute ``\\dim \\mathrm{HH}^n(X) = \\sum_{p+q=n} h^q(X, \\bigwedge^p T_X)``.
"""
function hochschild_dimension(P::PolyvectorParallelogram, n::Int)
  result = BigInt(0)
  for p in 0:P.dim
    q = n - p
    0 <= q <= P.dim || continue
    result += P[p, q]
  end
  result
end

# ─── Display as parallelogram ────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", P::PolyvectorParallelogram)
  d = P.dim

  # Width of largest number for alignment
  max_val = maximum(P.data)
  w = max(ndigits(max_val), 1)

  println(io, "Polyvector parallelogram (dim = $d):")
  println(io)

  # Row n goes from p+q = 0 to p+q = 2d
  # In row n: p ranges over max(0, n-d)..min(n, d)
  #           q = n - p
  for n in 0:(2 * d)
    # Indentation: (2d - n) spaces at top, 0 at middle, grows again
    # Actually the parallelogram shape: indent by (d - min(n, 2d-n)) etc.
    # The standard shape: each row n is indented by max(0, d - n) positions
    indent = max(0, d - n) * (w + 2)
    print(io, " " ^ indent)

    entries = String[]
    for p in 0:d
      q = n - p
      0 <= q <= d || continue
      push!(entries, lpad(string(P[p, q]), w))
    end

    println(io, join(entries, "  "))
  end
end

function Base.show(io::IO, P::PolyvectorParallelogram)
  d = P.dim
  print(io, "PolyvectorParallelogram(dim=$d)")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hochschild cohomology via HKR decomposition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hochschild_cohomology(X::PartialFlagVariety) -> PolyvectorParallelogram

Compute the Hochschild cohomology of ``X`` via the HKR decomposition:

``\\mathrm{HH}^n(X) = \\bigoplus_{p+q=n} H^q(X, \\bigwedge^p T_X)``

Returns a [`PolyvectorParallelogram`](@ref) encoding the full decomposition.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(2);

julia> P = hochschild_cohomology(X);

julia> P[0, 0]  # h⁰(𝒪) = 1
1

julia> P[2, 0]  # h⁰(∧²T) = h⁰(𝒪(3)) for ℙ²
10
```
"""
function hochschild_cohomology(X::PartialFlagVariety{MDT}) where {MDT}
  d = dimension(X)
  data = zeros(BigInt, d + 1, d + 1)

  T = tangent_bundle(X)

  for p in 0:d
    wedge_p = exterior_power(T, p)
    coh = dimensions(wedge_p)
    for q in 0:d
      data[p + 1, q + 1] = coh[q]
    end
  end

  PolyvectorParallelogram(data, d)
end
