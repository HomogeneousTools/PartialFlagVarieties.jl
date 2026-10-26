export volume, anticanonical_volume, invariants

_anticanonical_polarization(X::PartialFlagVariety) = anticanonical_bundle(X)
_anticanonical_polarization(Z::ZeroLocus) = tensor_product(
  anticanonical_bundle(ambient_variety(Z)), dual(det(defining_bundle(Z)))
)

"""
    volume(X::Union{PartialFlagVariety,ZeroLocus}) -> BigInt

Compute the anticanonical volume, the top self-intersection of the
anticanonical class. It is the leading coefficient of its Hilbert polynomial
multiplied by the factorial of the dimension, computed exactly.

For a zero locus, require [`is_strongly_fano`](@ref): its anticanonical class
must lift to an ample line bundle on the ambient variety. An `ArgumentError`
means that this sufficient condition was not certified, rather than a proof
that the zero locus is not Fano. Empty recorded zero loci have volume zero.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> volume(projective_space(3))
64

julia> volume(zero_locus(line_bundle(projective_space(4), 3)))
24
```
"""
function volume(X::Union{PartialFlagVariety,ZeroLocus})
  X isa ZeroLocus && !is_strongly_fano(X) &&
    throw(
      ArgumentError(
        "The anticanonical class must lift to an ample line bundle on the ambient variety."
      ),
    )
  polynomial = hilbert_polynomial(structure_sheaf(X), _anticanonical_polarization(X))
  d = dimension(X)
  coefficient = length(polynomial) > d ? polynomial[d + 1] : zero(Rational{BigInt})
  value = factorial(big(d)) * coefficient
  denominator(value) == 1 ||
    throw(ArgumentError("The anticanonical degree is not integral."))
  numerator(value)
end

"""
    anticanonical_volume(X::Union{PartialFlagVariety,ZeroLocus}) -> BigInt

Synonym for [`volume`](@ref).
"""
anticanonical_volume(X::Union{PartialFlagVariety,ZeroLocus}) = volume(X)

"""
    invariants(X::Union{PartialFlagVariety,ZeroLocus};
               hodge=true, max_tangent_power=min(1, dimension(X)))

Collect numerical invariants in a named tuple:

- `dimension`: the dimension of the variety or recorded zero locus;
- `volume`: [`volume`](@ref), or `nothing` when the ambient ampleness
  criterion for a zero locus is not satisfied;
- `chi_structure`: the holomorphic Euler characteristic of its structure sheaf;
- `chi_tangent`: Euler characteristics of exterior tangent powers, indexed
  by powers `1:max_tangent_power`;
- `chi_cotangent`: Euler characteristics of exterior cotangent powers,
  with entry `p+1` corresponding to power `p`, through `dimension(X) ÷ 2`;
- `hodge_numbers`: [`hodge_numbers`](@ref), or `nothing` if `hodge=false`.

Require `0 ≤ max_tangent_power ≤ dimension(X)`. Euler characteristics are
exact even when individual cohomology groups remain symbolic. The structure
sheaf Euler characteristic is holomorphic, including for a partial flag
variety, whose variety-valued `euler_characteristic(X)` is topological.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> I = invariants(projective_space(2); hodge=false);

julia> (I.dimension, I.volume, I.chi_structure)
(2, 9, 1)

julia> I.chi_tangent, I.chi_cotangent, I.hodge_numbers
(BigInt[8], BigInt[1, -1], nothing)
```
"""
function invariants(
  X::Union{PartialFlagVariety,ZeroLocus};
  hodge::Bool=true,
  max_tangent_power::Integer=min(1, dimension(X)),
)
  d = dimension(X)
  0 <= max_tangent_power <= d || throw(ArgumentError(
    "max_tangent_power must lie between zero and the dimension $d."
  ))
  # Euler characteristics depend only on the associated graded on a flag
  # variety; the ordinary tangent bundle also covers dimension zero.
  TX = tangent_bundle(X)
  Omega = dual(TX)
  chi_tangent = BigInt[
    euler_characteristic(exterior_power(TX, p)) for p in 1:max_tangent_power
  ]
  chi_cotangent = BigInt[
    euler_characteristic(exterior_power(Omega, p)) for p in 0:(d ÷ 2)
  ]
  fano_volume = X isa PartialFlagVariety || is_strongly_fano(X) ? volume(X) : nothing
  hodge_result = hodge ? hodge_numbers(X) : nothing
  (
    dimension=d,
    volume=fano_volume,
    chi_structure=euler_characteristic(structure_sheaf(X)),
    chi_tangent=chi_tangent,
    chi_cotangent=chi_cotangent,
    hodge_numbers=hodge_result,
  )
end
