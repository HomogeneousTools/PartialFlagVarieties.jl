export schur_functor

# Signed multiplicities are used only while evaluating Jacobi--Trudi.
function _virtual_tensor(a::Dict{IrrepLevi,BigInt}, b::Dict{IrrepLevi,BigInt})
  result = Dict{IrrepLevi,BigInt}()
  for (u, m) in a, (v, n) in b
    for (w, k) in _tensor_product_terms(u, v)
      result[w] = get(result, w, BigInt(0)) + m * n * k
    end
  end
  filter!(p -> !iszero(last(p)), result)
end

_virtual_bundle(E::CompletelyReducibleBundle) =
  Dict(u => BigInt(m) for (u, m) in _to_counts(E))

function _genuine_bundle(X::PartialFlagVariety, counts::Dict{IrrepLevi,BigInt})
  any(<(0), values(counts)) && error("Schur decomposition has negative multiplicities.")
  CompletelyReducibleBundle(X, IrrepLevi[u for (u, m) in counts for _ in 1:Int(m)])
end

# Unlike tensor powers, this also handles negative powers of a line bundle.
function _line_power(L::CompletelyReducibleBundle, k::Integer)
  rank(L) == 1 || throw(ArgumentError("Expected a line bundle."))
  CompletelyReducibleBundle(variety(L), k * p_dominant_weight(only(components(L))))
end

"""
    schur_functor(E::CompletelyReducibleBundle, λ)

Apply the Schur functor indexed by a weakly decreasing integer sequence `λ`.
Nonnegative partitions are padded with zeros to `rank(E)`; negative entries
are allowed only in a full length sequence and represent determinant twists.
A partition with more than `rank(E)` nonzero parts gives the zero bundle.

Works for direct sums as well as irreducible bundles, retaining their central
characters. Jacobi--Trudi uses the smaller of the row and column determinants;
determinant twists and duality reduce the polynomial degree first.
"""
function schur_functor(E::CompletelyReducibleBundle, λ::AbstractVector{<:Integer})
  issorted(λ; rev=true) || throw(ArgumentError("Schur weights must be weakly decreasing."))
  X = variety(E)
  r = Int(rank(E))
  a = Int.(λ)
  if any(<(0), a)
    length(a) == r ||
      throw(ArgumentError("Negative Schur weights require exactly rank(E) entries."))
  else
    while !isempty(a) && last(a) == 0
      pop!(a)
    end
    length(a) > r && return zero_bundle(X)
  end
  isempty(a) && return structure_sheaf(X)
  append!(a, zeros(Int, r - length(a)))
  shift = last(a)
  a .-= shift
  determinant_twist = _line_power(det(E), shift)
  all(iszero, a) && return determinant_twist

  # S_a(E) = det(E)^a₁ ⊗ S_(a₁-aᵣ,...,a₁-a₁)(E∨).
  complement = first(a) .- reverse(a)
  if sum(complement) < sum(a)
    return tensor_product(
      _line_power(det(E), shift + first(a)), schur_functor(dual(E), complement)
    )
  end
  filter!(!iszero, a)
  polynomial = if length(a) == 1
    symmetric_power(E, only(a))
  elseif first(a) == 1
    exterior_power(E, length(a))
  else
    _jacobi_trudi(E, a)
  end
  tensor_product(determinant_twist, polynomial)
end

function _jacobi_trudi(E::CompletelyReducibleBundle, a::Vector{Int})
  use_exterior = first(a) < length(a)
  λ = use_exterior ? [count(>=(j), a) for j in 1:first(a)] : a
  power = use_exterior ? exterior_power : symmetric_power
  n = length(λ)
  n < 8 * sizeof(UInt) ||
    throw(ArgumentError("Schur determinant is too large for subset evaluation."))
  powers = Dict(k => _virtual_bundle(power(E, k)) for k in 0:(first(λ) + n - 1))
  states = Dict(UInt(0) => _virtual_bundle(structure_sheaf(variety(E))))
  for i in 1:n
    next = Dict{UInt,Dict{IrrepLevi,BigInt}}()
    for (mask, value) in states, j in 1:n
      bit = UInt(1) << (j - 1)
      mask & bit == 0 || continue
      k = λ[i] - i + j
      k < 0 && continue
      sign = isodd(count_ones(mask >> j)) ? -1 : 1
      target = get!(() -> Dict{IrrepLevi,BigInt}(), next, mask | bit)
      for (u, m) in _virtual_tensor(value, powers[k])
        target[u] = get(target, u, BigInt(0)) + sign * m
      end
      filter!(p -> !iszero(last(p)), target)
    end
    states = next
  end
  _genuine_bundle(variety(E), states[(UInt(1) << n) - 1])
end

"""
    schur_functor(F::FilteredBundle, λ)

The induced Schur filtration. An auxiliary line-bundle weight records the
filtration degree during the Schur decomposition, so extensions are retained
even when Jacobi--Trudi cancels intermediate virtual summands.
"""
function schur_functor(F::FilteredBundle, λ::AbstractVector{<:Integer})
  X = variety(F)
  P = projective_space(1)
  A = product(X, P)
  tagged = zero_bundle(A)
  for (i, piece) in enumerate(graded_pieces(F))
    tagged = direct_sum(tagged, external_tensor_product(piece, line_bundle(P, i)))
  end
  result = schur_functor(tagged, λ)
  b = rank(dynkin_type(X))
  pieces = Dict{Int,Vector{IrrepLevi}}()
  for rep in components(result)
    c = collect(coefficients(p_dominant_weight(rep)))
    push!(get!(() -> IrrepLevi[], pieces, c[end]), IrrepLevi(marked_dynkin_type(X), c[1:b]))
  end
  FilteredBundle(
    X,
    [
      CompletelyReducibleBundle(X, pieces[degree]) for
      degree in sort!(collect(keys(pieces)))
    ],
  )
end
