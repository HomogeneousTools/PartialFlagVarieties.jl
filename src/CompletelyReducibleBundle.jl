# ═══════════════════════════════════════════════════════════════════════════════
#  CompletelyReducibleBundle — a completely reducible equivariant bundle on G/P
#
#  A completely reducible equivariant bundle on G/P is an equivariant bundle
#  whose underlying P-representation is completely reducible (semisimple).
#  It is encoded as a formal (virtual) sum of irreducible Levi representations,
#  together with a reference to the underlying partial flag variety.
# ═══════════════════════════════════════════════════════════════════════════════

export Bundle
export CompletelyReducibleBundle
export components, variety
export rank_bundle, tangent_bundle, cotangent_bundle
export structure_sheaf, O, zero_bundle, line_bundle, canonical_bundle, anticanonical_bundle
export det_bundle
export fano_index

# Names from Lie, StaticArrays, Combinatorics are available via the parent module.

# ═══════════════════════════════════════════════════════════════════════════════
#  Abstract Bundle type
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Bundle{MDT}

Abstract supertype for equivariant vector bundles on the partial flag variety
``G/P`` encoded by the marked Dynkin type `MDT`.

Concrete subtypes:
- [`CompletelyReducibleBundle`](@ref): semisimple equivariant bundles
- [`FilteredBundle`](@ref): bundles with a filtration by equivariant subbundles
"""
abstract type Bundle{MDT<:MarkedDynkinType} end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type definition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    CompletelyReducibleBundle{MDT}

A completely reducible equivariant vector bundle on the partial flag variety
``G/P`` encoded by the marked Dynkin type `MDT`.

Stored as a list of irreducible Levi representations (with multiplicity
encoded by repetition), together with the underlying [`PartialFlagVariety`](@ref).

# Fields
- `variety::PartialFlagVariety{MDT}`: the partial flag variety
- `components::Vector{IrrepLevi{MDT}}`: the irreducible summands

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> T = tangent_bundle(X);

julia> rank_bundle(T)
4
```
"""
struct CompletelyReducibleBundle{MDT<:MarkedDynkinType} <: Bundle{MDT}
  variety::PartialFlagVariety{MDT}
  components::Vector{IrrepLevi{MDT}}
end

"""
    CompletelyReducibleBundle(mdt::MarkedDynkinType, λ::WeightLatticeElem) -> CompletelyReducibleBundle

Construct the equivariant bundle on the partial flag variety determined by `mdt`
corresponding to the Levi representation of highest weight `λ`.
"""
function CompletelyReducibleBundle(mdt::MarkedDynkinType, λ::WeightLatticeElem)
  MDT = typeof(mdt)
  rep = IrrepLevi(MDT, λ)
  return CompletelyReducibleBundle(PartialFlagVariety{MDT}(), IrrepLevi{MDT}[rep])
end

function CompletelyReducibleBundle(X::PartialFlagVariety, λ::WeightLatticeElem)
  E = CompletelyReducibleBundle(marked_dynkin_type(X), λ)
  return CompletelyReducibleBundle(X, components(E))
end

# ─── Accessors ───────────────────────────────────────────────────────────────

"""
    variety(E::CompletelyReducibleBundle) -> PartialFlagVariety

Return the partial flag variety on which this bundle lives.
"""
variety(E::CompletelyReducibleBundle) = E.variety

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

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> rank_bundle(structure_sheaf(X))
1
```
"""
function rank_bundle(E::CompletelyReducibleBundle)
  return sum(fiber_dimension(c) for c in E.components; init=BigInt(0))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Standard bundles
# ═══════════════════════════════════════════════════════════════════════════════

"""
  O(X::PartialFlagVariety) -> CompletelyReducibleBundle
  structure_sheaf(X::PartialFlagVariety) -> CompletelyReducibleBundle

The trivial line bundle ``\\mathcal{O}_{G/P}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> rank_bundle(structure_sheaf(X))
1
```
"""
function O(X::PartialFlagVariety)
  return CompletelyReducibleBundle(X, WeightLatticeElem(dynkin_type(X)))
end

"""
    structure_sheaf(X::PartialFlagVariety) -> CompletelyReducibleBundle

The trivial line bundle ``\\mathcal{O}_{G/P}``. Alias for `O`.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> rank_bundle(structure_sheaf(X))
1
```
"""
structure_sheaf(X::PartialFlagVariety) = O(X)

"""
    zero_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The zero bundle on `X` (the empty direct sum).
"""
function zero_bundle(X::PartialFlagVariety)
  MDT = typeof(marked_dynkin_type(X))
  return CompletelyReducibleBundle{MDT}(X, IrrepLevi{MDT}[])
end

"""
  line_bundle(X::PartialFlagVariety, i::Integer) -> CompletelyReducibleBundle

The line bundle ``\\mathcal{O}(i)`` on `X`.

When the Picard rank of `X` is 1 (generalized Grassmannian), `i` is the
degree. When the Picard rank is greater than 1, this method raises an
error — use [`line_bundle(X, degrees)`](@ref) instead.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> L = line_bundle(X, 1);

julia> rank_bundle(L)
1
```
"""
function line_bundle(X::PartialFlagVariety, i::Integer)
  picard_rank(X) == 1 || throw(ArgumentError(
    "line_bundle(X, i::Integer) requires Picard rank 1, but X has Picard rank $(picard_rank(X)). " *
    "Use line_bundle(X, degrees::Vector{<:Integer}) instead."
  ))

  return line_bundle(X, [i])
end

"""
    line_bundle(X::PartialFlagVariety, degrees::Vector{<:Integer}) -> CompletelyReducibleBundle

The line bundle ``\\mathcal{O}(d_1, \\ldots, d_r)`` on `X`, where ``d_j`` is the
degree at the ``j``-th marked node.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{3}, (1, 3));

julia> L = line_bundle(X, [2, 1]);

julia> rank_bundle(L)
1
```
"""
function line_bundle(X::PartialFlagVariety, degrees::Vector{<:Integer})
  marked = marked_nodes(X)
  length(degrees) == length(marked) || throw(ArgumentError(
    string("Expected ", length(marked), " degrees (one per marked node), got ", length(degrees), ".")
  ))

  coords = zeros(Int, rank(X))
  for (j, m) in enumerate(marked)
    coords[m] = Int(degrees[j])
  end

  λ = WeightLatticeElem(dynkin_type(X), SVector{rank(X),Int}(coords))
  return CompletelyReducibleBundle(X, λ)
end

# ─── Type-level caches for tangent/cotangent rep lists ───────────────────────
# These avoid repeatedly computing tangent_weights + IrrepLevi decomposition
# for the same MDT.  The caches store only the IrrepLevi component vectors;
# fresh CompletelyReducibleBundle wrappers are created with the correct parent.

const _tangent_reps_cache = Dict{DataType,Vector}()
const _cotangent_reps_cache = Dict{DataType,Vector}()

function _tangent_reps(::Type{MDT}) where {MDT<:MarkedDynkinType}
  get!(_tangent_reps_cache, MDT) do
    tw = tangent_weights(MDT)
    IrrepLevi{MDT}[IrrepLevi(MDT, w) for w in tw]
  end::Vector{IrrepLevi{MDT}}
end

function _cotangent_reps(::Type{MDT}) where {MDT<:MarkedDynkinType}
  get!(_cotangent_reps_cache, MDT) do
    IrrepLevi{MDT}[dual(r) for r in _tangent_reps(MDT)]
  end::Vector{IrrepLevi{MDT}}
end

"""
    tangent_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The tangent bundle ``T_{G/P}``.

The tangent space at the identity coset decomposes as the direct sum of
root spaces for the positive nonparabolic roots. The tangent bundle
(as an equivariant bundle) is the semisimplification: for each maximal
nonparabolic root class, we get one irreducible Levi summand.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> T = tangent_bundle(X);

julia> n_components(T)
1

julia> rank_bundle(T)
4
```
"""
function tangent_bundle(X::PartialFlagVariety{MDT}) where {MDT<:MarkedDynkinType}
  return CompletelyReducibleBundle{MDT}(X, _tangent_reps(MDT))
end

"""
    cotangent_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The cotangent bundle ``\\Omega^1_{G/P} = T^*_{G/P}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> rank_bundle(cotangent_bundle(X)) == rank_bundle(tangent_bundle(X))
true
```
"""
function cotangent_bundle(X::PartialFlagVariety{MDT}) where {MDT<:MarkedDynkinType}
  return CompletelyReducibleBundle{MDT}(X, _cotangent_reps(MDT))
end

"""
    canonical_bundle(::Type{MDT}) -> CompletelyReducibleBundle
    canonical_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The canonical line bundle ``\\omega_{G/P} = K_{G/P}``.

This is computed directly from the formula

```math
K_{G/P} = -\\sum_{i \\in \\mathrm{marked}} a_i\\,\\omega_i,
\\qquad a_i = \\langle 2(\\rho_G - \\rho_P),\\,\\alpha_i^\\vee\\rangle,
```

without constructing the (co)tangent bundle.  See `anticanonical_degrees`.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = projective_space(1);

julia> K = canonical_bundle(X);

julia> rank_bundle(K)
1

julia> dimensions(K)[1]  # H¹(ℙ¹, 𝒪(-2)) = 1
1

julia> K2 = canonical_bundle(projective_space(4));

julia> dimensions(K2)[4]  # H⁴(ℙ⁴, 𝒪(-5)) = 1
1
```
"""
function canonical_bundle(X::PartialFlagVariety{MDT}) where {MDT<:MarkedDynkinType}
  degs = anticanonical_degrees(MDT)
  return line_bundle(X, Vector{Int}(-degs))
end

"""Type-level dispatch of [`canonical_bundle`](@ref): returns a zero-variety
bundle placeholder for use in type-level computations.

Same as `canonical_bundle(X)` where `X = PartialFlagVariety{MDT}()`."""
function canonical_bundle(::Type{MDT}) where {MDT<:MarkedDynkinType}
  X = PartialFlagVariety{MDT}()
  canonical_bundle(X)
end

"""
    anticanonical_bundle(::Type{MDT}) -> CompletelyReducibleBundle
    anticanonical_bundle(X::PartialFlagVariety) -> CompletelyReducibleBundle

The anticanonical line bundle ``\\omega_{G/P}^{-1} = -K_{G/P}``.

This is computed directly from the formula

```math
-K_{G/P} = \\sum_{i \\in \\mathrm{marked}} a_i\\,\\omega_i,
\\qquad a_i = \\langle 2(\\rho_G - \\rho_P),\\,\\alpha_i^\\vee\\rangle.
```

See `anticanonical_degrees`.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = projective_space(1);

julia> L = anticanonical_bundle(X);

julia> rank_bundle(L)
1

julia> dimensions(L)[0]  # H⁰(ℙ¹, 𝒪(2)) = 3
3

julia> L2 = anticanonical_bundle(projective_space(4));

julia> dimensions(L2)[0]  # H⁰(ℙ⁴, 𝒪(5)) = 126
126
```
"""
function anticanonical_bundle(X::PartialFlagVariety{MDT}) where {MDT<:MarkedDynkinType}
  degs = anticanonical_degrees(MDT)
  return line_bundle(X, Vector{Int}(degs))
end

"""Type-level dispatch of [`anticanonical_bundle`](@ref): returns a zero-variety
bundle placeholder for use in type-level computations.

Same as `anticanonical_bundle(X)` where `X = PartialFlagVariety{MDT}()`."""
function anticanonical_bundle(::Type{MDT}) where {MDT<:MarkedDynkinType}
  X = PartialFlagVariety{MDT}()
  anticanonical_bundle(X)
end

"""
    anticanonical_degrees(X::PartialFlagVariety) -> SVector

Instance-level dispatch of [`anticanonical_degrees(::Type{MDT})`](@ref).
"""
anticanonical_degrees(X::PartialFlagVariety{MDT}) where {MDT} = anticanonical_degrees(MDT)

"""
    fano_index(X::PartialFlagVariety) -> Int

The Fano index of the partial flag variety ``G/P``, defined as the unique
positive integer ``r`` such that

```math
-K_{G/P} = r\\,\\omega_m
```

where ``\\omega_m`` is the ample generator of ``\\mathrm{Pic}(G/P) \\cong \\mathbb{Z}``.

All partial flag varieties are Fano (the anticanonical bundle is ample), so
the Fano index is always a positive integer.

Throws an `ArgumentError` when `picard_rank(X) > 1`; use
`anticanonical_degrees` for the multi-degree description in that case.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> fano_index(projective_space(4))
5

julia> fano_index(Gr(2, 5))
5

julia> fano_index(quadric(4))
4
```
"""
function fano_index(X::PartialFlagVariety{MDT}) where {MDT<:MarkedDynkinType}
  Marked = marked_nodes(MDT)
  length(Marked) == 1 || throw(ArgumentError(
    "fano_index is only defined for Picard-rank-1 varieties; " *
    "use anticanonical_degrees for the multi-degree description.")
  )
  anticanonical_degrees(MDT)[1]
end

function _trivial_semisimple_weight(X::PartialFlagVariety)
  MDT = typeof(marked_dynkin_type(X))
  LT = levi_type(MDT)
  return WeightLatticeElem(LT === nothing ? TypeA{1} : LT)
end

function _det_bundle_irrep(X::PartialFlagVariety, rep::IrrepLevi)
  MDT = typeof(marked_dynkin_type(X))
  d = Int(fiber_dimension(rep))
  return IrrepLevi(MDT, d * rep.central, _trivial_semisimple_weight(X))
end

"""
    det_bundle(E::CompletelyReducibleBundle) -> CompletelyReducibleBundle

The determinant line bundle ``\\det(E) = \\bigwedge^{\\mathrm{rk}(E)} E``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> rank_bundle(det_bundle(tangent_bundle(X)))
1
```
"""
function det_bundle(E::CompletelyReducibleBundle)
  X = variety(E)
  λ = WeightLatticeElem(dynkin_type(X))
  for c in components(E)
    λ += p_dominant_weight(_det_bundle_irrep(X, c))
  end
  return CompletelyReducibleBundle(X, λ)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Monoidal operations on bundles
# ═══════════════════════════════════════════════════════════════════════════════

"""
    dual(E::CompletelyReducibleBundle) -> CompletelyReducibleBundle

The dual bundle ``E^*``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = Gr(2, 4);

julia> E = tangent_bundle(X);

julia> rank_bundle(dual(E)) == rank_bundle(E)
true
```
"""
function dual(E::CompletelyReducibleBundle{MDT}) where {MDT}
  return CompletelyReducibleBundle{MDT}(E.variety, [dual(c) for c in E.components])
end

"""
    tensor_product(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)

The tensor product ``E \\otimes F``.

Uses bilinearity: ``(\\bigoplus_i V_i) \\otimes (\\bigoplus_j W_j) = \\bigoplus_{i,j} V_i \\otimes W_j``

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = tangent_bundle(X);

julia> S = structure_sheaf(X);

julia> rank_bundle(tensor_product(E, S)) == rank_bundle(E)
true
```
"""
function tensor_product(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  X = variety(E)
  Y = variety(F)
  MDT = typeof(marked_dynkin_type(X))
  typeof(marked_dynkin_type(Y)) === MDT || throw(ArgumentError(
    "tensor_product requires bundles on the same partial flag variety type."
  ))

  # Deduplicate components to avoid redundant tensor product calls.
  # E.g., O(1)⁶ has 6 identical components → compute ⊗ once, replicate.
  e_counts = Dict{IrrepLevi{MDT},Int}()
  for a in E.components
    e_counts[a] = get(e_counts, a, 0) + 1
  end
  f_counts = Dict{IrrepLevi{MDT},Int}()
  for b in F.components
    f_counts[b] = get(f_counts, b, 0) + 1
  end

  result = IrrepLevi{MDT}[]
  for (a, ma) in e_counts
    for (b, mb) in f_counts
      tp = tensor_product(a, b)
      total = ma * mb
      for _ in 1:total
        append!(result, tp)
      end
    end
  end
  return CompletelyReducibleBundle{MDT}(X, result)
end

"""
    exterior_power(E::CompletelyReducibleBundle, k::Int) -> CompletelyReducibleBundle

The k-th exterior power ``\\bigwedge^k E``.

For a direct sum ``E = \\bigoplus_i V_i``, we use:
``\\bigwedge^k E = \\bigoplus_{|\\alpha| = k} \\bigotimes_i \\bigwedge^{\\alpha_i} V_i``

where ``\\alpha`` runs over compositions (multiexponents) of ``k`` in
``n = n_{\\rm components}`` parts.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = tangent_bundle(X);

julia> rank_bundle(exterior_power(E, 0))
1

julia> rank_bundle(exterior_power(E, 1)) == rank_bundle(E)
true
```
"""
function exterior_power(E::CompletelyReducibleBundle{MDT}, k::Integer) where {MDT}
  k = Int(k)
  n = n_components(E)
  k < 0 && return CompletelyReducibleBundle{MDT}(E.variety, IrrepLevi{MDT}[])
  k == 0 && return structure_sheaf(E.variety)
  k == 1 && return E

  ranks = [Int(fiber_dimension(c)) for c in E.components]

  result = IrrepLevi{MDT}[]

  for α in multiexponents(n, k)
    any(α[i] > ranks[i] for i in 1:n) && continue

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

  return CompletelyReducibleBundle{MDT}(E.variety, result)
end

"""
    symmetric_power(E::CompletelyReducibleBundle, k::Int) -> CompletelyReducibleBundle

The k-th symmetric power ``\\mathrm{Sym}^k E``.

For a direct sum ``E = \\bigoplus_i V_i``, we use:
``\\mathrm{Sym}^k E = \\bigoplus_{|\\alpha| = k} \\bigotimes_i \\mathrm{Sym}^{\\alpha_i} V_i``

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = tangent_bundle(X);

julia> rank_bundle(symmetric_power(E, 0))
1

julia> rank_bundle(symmetric_power(E, 1)) == rank_bundle(E)
true
```
"""
function symmetric_power(E::CompletelyReducibleBundle{MDT}, k::Integer) where {MDT}
  k = Int(k)
  n = n_components(E)
  k < 0 && return CompletelyReducibleBundle{MDT}(E.variety, IrrepLevi{MDT}[])
  k == 0 && return structure_sheaf(E.variety)
  k == 1 && return E

  result = IrrepLevi{MDT}[]

  for α in multiexponents(n, k)
    factors = Vector{Vector{IrrepLevi{MDT}}}()
    skip = false
    for i in 1:n
      sym_i = symmetric_power(E.components[i], α[i])
      # TODO should be iszero?
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

  return CompletelyReducibleBundle{MDT}(E.variety, result)
end

"""
    direct_sum(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)

The direct sum ``E \\oplus F``.
"""
function direct_sum(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle)
  X = variety(E)
  Y = variety(F)
  MDT = typeof(marked_dynkin_type(X))
  typeof(marked_dynkin_type(Y)) === MDT || throw(ArgumentError(
    "direct_sum requires bundles on the same partial flag variety type."
  ))
  return CompletelyReducibleBundle{MDT}(X, vcat(E.components, F.components))
end

# ─── Twist ───────────────────────────────────────────────────────────────────

"""
    twist(E::CompletelyReducibleBundle, i::Int, k::Int=1) -> CompletelyReducibleBundle

Twist ``E`` by ``\\mathcal{O}(k)`` at the `i`-th marked node:
``E(k) = E \\otimes \\mathcal{L}(k \\omega_{m_i})``.

# Examples
```jldoctest
julia> using PartialFlagVarieties, Lie

julia> X = partial_flag_variety(TypeA{4}, (1,));

julia> E = structure_sheaf(X);

julia> rank_bundle(twist(E, 1, 3))
1
```
"""
function twist(E::CompletelyReducibleBundle{MDT}, i::Integer, k::Integer=1) where {MDT<:MarkedDynkinType{DT,Marked}} where {DT,Marked}
  i, k = Int(i), Int(k)
  1 <= i <= length(Marked) || throw(ArgumentError(
    "Index $i out of range. MDT has $(length(Marked)) marked node(s)."
  ))

  m = Marked[i]
  ω = fundamental_weight(DT, m)
  λ = k * ω
  twist_rep = IrrepLevi(MDT, λ)

  result = IrrepLevi{MDT}[]
  for c in E.components
    append!(result, tensor_product(c, twist_rep))
  end
  return CompletelyReducibleBundle{MDT}(E.variety, result)
end

# ─── Arithmetic operators ───────────────────────────────────────────────────

Base.:*(E::CompletelyReducibleBundle, F::CompletelyReducibleBundle) = tensor_product(E, F)

"""
    n * E -> CompletelyReducibleBundle

The `n`-fold direct sum ``E \\oplus \\cdots \\oplus E`` (n ≥ 1).
Returns the zero bundle when `n == 0` and throws for `n < 0`.
"""
function Base.:*(n::Integer, E::CompletelyReducibleBundle{MDT}) where {MDT}
  n < 0 && throw(ArgumentError("Cannot multiply a bundle by a negative integer ($n)"))
  n == 0 && return zero_bundle(E.variety)
  return CompletelyReducibleBundle{MDT}(E.variety, repeat(E.components, n))
end
Base.:*(E::CompletelyReducibleBundle, n::Integer) = n * E

const ⊕ = direct_sum
const ⊗ = tensor_product

export ⊕, ⊗

# ─── iszero ──────────────────────────────────────────────────────────

"""
    iszero(E::CompletelyReducibleBundle) -> Bool

Return `true` if `E` is the zero bundle (has no summands).
"""
Base.iszero(E::CompletelyReducibleBundle) = isempty(E.components)

# ─── Display ─────────────────────────────────────────────────────────────────

function Base.show(io::IO, E::CompletelyReducibleBundle)
  if isempty(E.components)
    print(io, "0")
  elseif length(E.components) == 1
    print(io, "E", E.components[1])
  else
    parts = ["E" * sprint(show, c) for c in E.components]
    print(io, join(parts, " ⊕ "))
  end
end
