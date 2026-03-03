# ═══════════════════════════════════════════════════════════════════════════════
#  CompletelyReducibleBundle — a completely reducible equivariant bundle on G/P
#
#  A completely reducible equivariant bundle on G/P is an equivariant bundle
#  whose underlying P-representation is completely reducible (semisimple).
#  It is encoded as a formal (virtual) sum of irreducible Levi representations.
# ═══════════════════════════════════════════════════════════════════════════════

export CompletelyReducibleBundle
export components
export rank_bundle, tangent_bundle, cotangent_bundle
export structure_sheaf, line_bundle, canonical_bundle, anticanonical_bundle
export det_bundle

# Names from Lie, StaticArrays, Combinatorics are available via the parent module.

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    CompletelyReducibleBundle{MDT}

A completely reducible equivariant vector bundle on the partial flag variety
``G/P`` encoded by the marked Dynkin type `MDT`.

Stored as a list of irreducible Levi representations (with multiplicity
encoded by repetition).

# Fields
- `components::Vector{IrrepLevi{MDT}}`: the irreducible summands

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> T = tangent_bundle(MDT);

julia> rank_bundle(T)
4
```
"""
struct CompletelyReducibleBundle{MDT<:MarkedDynkinType}
  components::Vector{IrrepLevi{MDT}}
end

# ─── Accessors ───────────────────────────────────────────────────────────────

"""
    components(E::CompletelyReducibleBundle) -> Vector{IrrepLevi}

Return the irreducible summands.
"""
components(E::CompletelyReducibleBundle) = E.components

"""
    n_components(E::CompletelyReducibleBundle) -> Int

Return the number of irreducible summands.
"""
n_components(E::CompletelyReducibleBundle) = length(E.components)

"""
    rank_bundle(E::CompletelyReducibleBundle) -> BigInt

Return the total rank (fiber dimension) of the bundle.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (1,)};

julia> rank_bundle(structure_sheaf(MDT))
1
```
"""
function rank_bundle(E::CompletelyReducibleBundle{MDT}) where {MDT}
  return sum(fiber_dimension(c) for c in E.components; init=BigInt(0))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Standard bundles
# ═══════════════════════════════════════════════════════════════════════════════

"""
    structure_sheaf(::Type{MDT}) -> CompletelyReducibleBundle{MDT}

The trivial line bundle ``\\mathcal{O}_{G/P}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> rank_bundle(structure_sheaf(MDT))
1
```
"""
function structure_sheaf(::Type{MDT}) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  R = rank(DT)
  LT = levi_type(MDT)

  zero_central = zeros(Rational{Int}, length(Marked))
  if LT === nothing
    # Full flag: use a dummy 0-dim semisimple part
    zero_ss = WeightLatticeElem(TypeA{1}, [0])
  else
    LR = rank(LT)
    zero_ss = WeightLatticeElem(LT, zeros(Int, LR))
  end

  triv = IrrepLevi{MDT}(zero_central, zero_ss)
  return CompletelyReducibleBundle{MDT}([triv])
end

"""
    line_bundle(::Type{MDT}, i::Int) -> CompletelyReducibleBundle{MDT}

The line bundle ``\\mathcal{O}(1)`` corresponding to the `i`-th marked node.
This is ``\\mathcal{L}(\\omega_{m_i})`` where ``m_i`` is the `i`-th marked node.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (2,)};

julia> L = line_bundle(MDT, 1);

julia> rank_bundle(L)
1
```
"""
function line_bundle(::Type{MDT}, i::Int=1) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  1 <= i <= length(Marked) || throw(ArgumentError(
    "Index $i out of range. MDT has $(length(Marked)) marked node(s)."
  ))
  m = Marked[i]
  ω = fundamental_weight(DT, m)
  rep = IrrepLevi(MDT, ω)
  return CompletelyReducibleBundle{MDT}([rep])
end

"""
    tangent_bundle(::Type{MDT}) -> CompletelyReducibleBundle{MDT}

The tangent bundle ``T_{G/P}``.

The tangent space at the identity coset decomposes as the direct sum of
root spaces for the positive nonparabolic roots. The tangent bundle
(as an equivariant bundle) is the semisimplification: for each maximal
nonparabolic root class, we get one irreducible Levi summand.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> T = tangent_bundle(MDT);

julia> n_components(T)
2

julia> rank_bundle(T)
4
```
"""
function tangent_bundle(::Type{MDT}) where {MDT<:MarkedDynkinType}
  tw = tangent_weights(MDT)
  reps = [IrrepLevi(MDT, w) for w in tw]
  return CompletelyReducibleBundle{MDT}(reps)
end

"""
    cotangent_bundle(::Type{MDT}) -> CompletelyReducibleBundle{MDT}

The cotangent bundle ``\\Omega^1_{G/P} = T^*_{G/P}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> rank_bundle(cotangent_bundle(MDT)) == rank_bundle(tangent_bundle(MDT))
true
```
"""
function cotangent_bundle(::Type{MDT}) where {MDT<:MarkedDynkinType}
  return dual(tangent_bundle(MDT))
end

"""
    canonical_bundle(::Type{MDT}) -> CompletelyReducibleBundle{MDT}

The canonical line bundle ``\\omega_{G/P} = \\det(\\Omega^1_{G/P})``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (1,)};

julia> rank_bundle(canonical_bundle(MDT))
1
```
"""
function canonical_bundle(::Type{MDT}) where {MDT<:MarkedDynkinType}
  return det_bundle(cotangent_bundle(MDT))
end

"""
    anticanonical_bundle(::Type{MDT}) -> CompletelyReducibleBundle{MDT}

The anticanonical line bundle ``\\omega_{G/P}^{-1}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (1,)};

julia> rank_bundle(anticanonical_bundle(MDT))
1
```
"""
function anticanonical_bundle(::Type{MDT}) where {MDT<:MarkedDynkinType}
  return dual(canonical_bundle(MDT))
end

"""
    det_bundle(E::CompletelyReducibleBundle{MDT}) -> CompletelyReducibleBundle{MDT}

The determinant line bundle ``\\det(E) = \\bigwedge^{\\mathrm{rk}(E)} E``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> rank_bundle(det_bundle(tangent_bundle(MDT)))
1
```
"""
function det_bundle(E::CompletelyReducibleBundle{MDT}) where {MDT}
  total_central = zeros(Rational{Int}, central_rank(MDT))
  for c in E.components
    # Each component contributes its central part scaled by its semisimple dimension
    d = fiber_dimension(c)
    total_central .+= d * c.central
  end

  LT = levi_type(MDT)

  # Determinant representation of the semisimple part
  # For each irreducible V_λ, det(V_λ) contributes one copy of the
  # sum of weights = dim(V_λ) * (sum of weights / dim)
  # But the determinant of a direct sum is the tensor product of determinants
  if LT === nothing
    zero_ss = WeightLatticeElem(TypeA{1}, [0])
  else
    LR = rank(LT)
    # The determinant of V_λ is the character that maps the maximal torus
    # to the wedge product. For a representation with highest weight λ,
    # det(V_λ) has weight = sum of all weights of V_λ.
    # A simpler approach: det of a direct sum ⊕ V_i = ⊗ det(V_i)
    # For the central part, det(V_λ) = dim(V_λ) * central_weight
    # For the semisimple part, det(V_λ) = 0 (it's a character of the center)
    zero_ss = WeightLatticeElem(LT, zeros(Int, LR))
  end

  return CompletelyReducibleBundle{MDT}([IrrepLevi{MDT}(total_central, zero_ss)])
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Monoidal operations on bundles
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dual(E::CompletelyReducibleBundle{MDT}) -> CompletelyReducibleBundle{MDT}

The dual bundle ``E^*``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{3}, (2,)};

julia> E = tangent_bundle(MDT);

julia> rank_bundle(dual(E)) == rank_bundle(E)
true
```
"""
function dual(E::CompletelyReducibleBundle{MDT}) where {MDT}
  return CompletelyReducibleBundle{MDT}([dual(c) for c in E.components])
end

"""
    tensor_product(E::CompletelyReducibleBundle{MDT}, F::CompletelyReducibleBundle{MDT})

The tensor product ``E \\otimes F``.

Uses bilinearity: ``(\\bigoplus_i V_i) \\otimes (\\bigoplus_j W_j) = \\bigoplus_{i,j} V_i \\otimes W_j``

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (1,)};

julia> E = tangent_bundle(MDT);

julia> S = structure_sheaf(MDT);

julia> rank_bundle(tensor_product(E, S)) == rank_bundle(E)
true
```
"""
function tensor_product(E::CompletelyReducibleBundle{MDT}, F::CompletelyReducibleBundle{MDT}) where {MDT}
  result = IrrepLevi{MDT}[]
  for a in E.components
    for b in F.components
      append!(result, tensor_product(a, b))
    end
  end
  return CompletelyReducibleBundle{MDT}(result)
end

"""
    direct_sum(E::CompletelyReducibleBundle{MDT}, F::CompletelyReducibleBundle{MDT})

The direct sum ``E \\oplus F``.
"""
function direct_sum(E::CompletelyReducibleBundle{MDT}, F::CompletelyReducibleBundle{MDT}) where {MDT}
  return CompletelyReducibleBundle{MDT}(vcat(E.components, F.components))
end

"""
    exterior_power(E::CompletelyReducibleBundle{MDT}, k::Int) -> CompletelyReducibleBundle{MDT}

The k-th exterior power ``\\bigwedge^k E``.

For a direct sum ``E = \\bigoplus_i V_i``, we use:
``\\bigwedge^k E = \\bigoplus_{|\\alpha| = k} \\bigotimes_i \\bigwedge^{\\alpha_i} V_i``

where ``\\alpha`` runs over compositions (multiexponents) of ``k`` in
``n = n_{\\rm components}`` parts.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (1,)};

julia> E = tangent_bundle(MDT);

julia> rank_bundle(exterior_power(E, 0))
1

julia> rank_bundle(exterior_power(E, 1)) == rank_bundle(E)
true
```
"""
function exterior_power(E::CompletelyReducibleBundle{MDT}, k::Int) where {MDT}
  n = n_components(E)
  k < 0 && return CompletelyReducibleBundle{MDT}(IrrepLevi{MDT}[])
  k == 0 && return structure_sheaf(MDT)
  k == 1 && return E

  # Compute fiber ranks of each irreducible summand
  ranks = [Int(fiber_dimension(c)) for c in E.components]

  result = IrrepLevi{MDT}[]

  # Iterate over multiexponents: compositions of k into n parts
  # where each part αᵢ ≤ rank(Vᵢ)
  for α in multiexponents(n, k)
    # Check feasibility: αᵢ ≤ rank of i-th summand
    any(α[i] > ranks[i] for i in 1:n) && continue

    # Build tensor product of ⋀^{αᵢ}(Vᵢ)
    factors = Vector{Vector{IrrepLevi{MDT}}}()
    skip = false
    for i in 1:n
      wedge_i = exterior_power(E.components[i], α[i])
      if isempty(wedge_i)
        skip = true
        break
      end
      push!(factors, wedge_i)
    end
    skip && continue

    # Tensor all factors together
    current = factors[1]
    for i in 2:length(factors)
      next = IrrepLevi{MDT}[]
      for a in current
        for b in factors[i]
          append!(next, tensor_product(a, b))
        end
      end
      current = next
    end
    append!(result, current)
  end

  return CompletelyReducibleBundle{MDT}(result)
end

"""
    symmetric_power(E::CompletelyReducibleBundle{MDT}, k::Int) -> CompletelyReducibleBundle{MDT}

The k-th symmetric power ``\\mathrm{Sym}^k E``.

For a direct sum ``E = \\bigoplus_i V_i``, we use:
``\\mathrm{Sym}^k E = \\bigoplus_{|\\alpha| = k} \\bigotimes_i \\mathrm{Sym}^{\\alpha_i} V_i``

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (1,)};

julia> E = tangent_bundle(MDT);

julia> rank_bundle(symmetric_power(E, 0))
1

julia> rank_bundle(symmetric_power(E, 1)) == rank_bundle(E)
true
```
"""
function symmetric_power(E::CompletelyReducibleBundle{MDT}, k::Int) where {MDT}
  n = n_components(E)
  k < 0 && return CompletelyReducibleBundle{MDT}(IrrepLevi{MDT}[])
  k == 0 && return structure_sheaf(MDT)
  k == 1 && return E

  result = IrrepLevi{MDT}[]

  for α in multiexponents(n, k)
    factors = Vector{Vector{IrrepLevi{MDT}}}()
    skip = false
    for i in 1:n
      sym_i = symmetric_power(E.components[i], α[i])
      if isempty(sym_i)
        skip = true
        break
      end
      push!(factors, sym_i)
    end
    skip && continue

    current = factors[1]
    for i in 2:length(factors)
      next = IrrepLevi{MDT}[]
      for a in current
        for b in factors[i]
          append!(next, tensor_product(a, b))
        end
      end
      current = next
    end
    append!(result, current)
  end

  return CompletelyReducibleBundle{MDT}(result)
end

# ─── Twist ───────────────────────────────────────────────────────────────────

"""
    twist(E::CompletelyReducibleBundle{MDT}, i::Int, k::Int=1) -> CompletelyReducibleBundle{MDT}

Twist ``E`` by ``\\mathcal{O}(k)`` at the `i`-th marked node:
``E(k) = E \\otimes \\mathcal{L}(k \\omega_{m_i})``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> MDT = MarkedDynkinType{TypeA{4}, (1,)};

julia> E = structure_sheaf(MDT);

julia> rank_bundle(twist(E, 1, 3))
1
```
"""
function twist(E::CompletelyReducibleBundle{MDT}, i::Int, k::Int=1) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  1 <= i <= length(Marked) || throw(ArgumentError(
    "Index $i out of range. MDT has $(length(Marked)) marked node(s)."
  ))

  # Create line bundle L(k·ω_{m_i})
  m = Marked[i]
  ω = fundamental_weight(DT, m)
  λ = k * ω
  twist_rep = IrrepLevi(MDT, λ)

  # Twist each component
  result = IrrepLevi{MDT}[]
  for c in E.components
    append!(result, tensor_product(c, twist_rep))
  end
  return CompletelyReducibleBundle{MDT}(result)
end

# ─── Arithmetic operators ───────────────────────────────────────────────────

Base.:*(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle) = tensor_product(E, F)

const ⊕ = direct_sum
const ⊗ = tensor_product

export ⊕, ⊗

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, E::CompletelyReducibleBundle{MDT}) where {MDT}
  if isempty(E.components)
    print(io, "0")
  elseif length(E.components) == 1
    print(io, "E", E.components[1])
  else
    parts = ["E" * sprint(show, c) for c in E.components]
    print(io, join(parts, " ⊕ "))
  end
end
