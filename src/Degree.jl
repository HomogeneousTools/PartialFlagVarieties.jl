"""
    degree(X::Union{PartialFlagVariety,ZeroLocus}) -> BigInt
    degree(X::Union{PartialFlagVariety,ZeroLocus}, L::CompletelyReducibleBundle) -> BigInt

Compute the degree of `X` with respect to a line bundle `L`, namely the top
self-intersection ``\\int_X c_1(L)^d``, where ``d = \\dim X``.
When `L` is very ample, this is the projective degree in the embedding it defines.

The default polarization is ``\\mathcal{O}(1, \\ldots, 1)``: degree one at every
marked node of `X`, or of its ambient variety for a zero locus.
On a Grassmannian this is the Plücker polarization.
It need not be the anticanonical polarization, or a primitive generator of the
Picard group of a zero locus.

An explicit `L` must be a completely reducible line bundle on `X`, or on its
ambient variety for a zero locus, as for [`hilbert_polynomial`](@ref).
No ampleness assumption is needed for the top self-intersection.
For a nonample line bundle the result can be zero or negative, and need not
equal its volume.

The result is computed exactly as ``d!`` times the coefficient of ``t^d`` in
``t \\mapsto \\chi(X, L^{\\otimes t})``, using the symbolic Hilbert polynomial.
The coefficient is zero when the polynomial has degree smaller than ``d``.
Dimension-zero loci return their scheme-theoretic length; empty recorded zero
loci return zero. As for the rest of the zero-locus API, regularity is assumed.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> degree(Gr(2, 5))
5

julia> X = projective_space(1) * projective_space(2);

julia> degree(X)  # Segre polarization O(1,1)
3

julia> degree(X, line_bundle(X, [2, 3]))
54

julia> G = Gr(3, 7);

julia> A = exterior_power(dual(universal_subbundle(G)), 2);

julia> degree(zero_locus(A ⊕ A))  # the catalecticant sixfold
61
```
"""
function degree(X::Union{PartialFlagVariety,ZeroLocus})
  ambient = X isa ZeroLocus ? ambient_variety(X) : X
  degree(X, line_bundle(ambient, ones(Int, picard_rank(ambient))))
end

function degree(X::Union{PartialFlagVariety,ZeroLocus}, L::CompletelyReducibleBundle)
  d = dimension(X)
  polynomial = hilbert_polynomial(structure_sheaf(X), L)
  value = factorial(big(d)) * polynomial[d]
  denominator(value) == 1 || error("The top self-intersection is not integral.")
  numerator(value)
end
