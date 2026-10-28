export Polynomial,
  coefficients, coefficient, constant_coefficient, leading_coefficient, degree

"""
    Polynomial(coefficients::AbstractVector)
    Polynomial(constant::Union{Integer,Rational})

A polynomial in `t` with exact `Rational{BigInt}` coefficients.
The vector constructor takes coefficients in increasing degree order and
discards trailing zeros; an empty vector defines the zero polynomial.
Only integer and rational coefficients are accepted.

Polynomials support addition, subtraction, multiplication, nonnegative integer
powers, and multiplication or division by rational scalars.
Evaluate a polynomial `p` at `t` with `p(t)`.
Use [`coefficient`](@ref), [`constant_coefficient`](@ref), and
[`leading_coefficient`](@ref) to access individual coefficients, or
[`coefficients(::Polynomial)`](@ref) to obtain a coefficient vector.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> p = Polynomial([1, 2, 1])
t^2 + 2*t + 1

julia> p(3)
16//1

julia> coefficient(p, 1)
2//1
```
"""
struct Polynomial
  coeffs::Vector{Rational{BigInt}}

  function Polynomial(coefficients::AbstractVector{<:Union{Integer,Rational}})
    coeffs = Rational{BigInt}.(coefficients)
    all(isfinite, coeffs) || throw(ArgumentError("Polynomial coefficients must be finite."))
    last_nonzero = findlast(!iszero, coeffs)
    if last_nonzero === nothing
      coeffs = Rational{BigInt}[0]
    else
      resize!(coeffs, last_nonzero)
    end
    new(coeffs)
  end
end

Polynomial(constant::Union{Integer,Rational}) = Polynomial([constant])
Polynomial() = Polynomial(0)

"""Return a copy of the coefficients of `p` in increasing degree order."""
coefficients(p::Polynomial) = copy(p.coeffs)

"""Return the coefficient of `t^k`, or zero above the degree; require `k ≥ 0`."""
function coefficient(p::Polynomial, k::Integer)
  k >= 0 || throw(DomainError(k, "The coefficient degree must be nonnegative."))
  k < length(p.coeffs) ? p.coeffs[Int(k) + 1] : zero(Rational{BigInt})
end

"""Return the coefficient of `t^0`."""
constant_coefficient(p::Polynomial) = first(p.coeffs)

"""Return the coefficient of highest degree, or zero for the zero polynomial."""
leading_coefficient(p::Polynomial) = last(p.coeffs)

"""Return the degree of `p`, with degree `-1` for the zero polynomial."""
degree(p::Polynomial) = iszero(p) ? -1 : length(p.coeffs) - 1

(p::Polynomial)(t::Number) = evalpoly(t, p.coeffs)

Base.iszero(p::Polynomial) = length(p.coeffs) == 1 && iszero(first(p.coeffs))
Base.zero(::Type{Polynomial}) = Polynomial(0)
Base.one(::Type{Polynomial}) = Polynomial(1)
Base.zero(::Polynomial) = zero(Polynomial)
Base.one(::Polynomial) = one(Polynomial)
Base.copy(p::Polynomial) = Polynomial(p.coeffs)
Base.:(==)(p::Polynomial, q::Polynomial) = p.coeffs == q.coeffs
Base.hash(p::Polynomial, h::UInt) = hash(p.coeffs, h)

function Base.:+(p::Polynomial, q::Polynomial)
  Polynomial([
    coefficient(p, k) + coefficient(q, k) for k in 0:max(degree(p), degree(q), 0)
  ])
end

Base.:-(p::Polynomial) = Polynomial(-p.coeffs)
Base.:-(p::Polynomial, q::Polynomial) = p + (-q)

function Base.:*(p::Polynomial, q::Polynomial)
  (iszero(p) || iszero(q)) && return zero(Polynomial)
  result = zeros(Rational{BigInt}, degree(p) + degree(q) + 1)
  for i in eachindex(p.coeffs), j in eachindex(q.coeffs)
    result[i + j - 1] += p.coeffs[i] * q.coeffs[j]
  end
  Polynomial(result)
end

Base.:+(p::Polynomial, c::Union{Integer,Rational}) = p + Polynomial(c)
Base.:+(c::Union{Integer,Rational}, p::Polynomial) = p + c
Base.:-(p::Polynomial, c::Union{Integer,Rational}) = p - Polynomial(c)
Base.:-(c::Union{Integer,Rational}, p::Polynomial) = Polynomial(c) - p
Base.:*(c::Union{Integer,Rational}, p::Polynomial) = Polynomial(c .* p.coeffs)
Base.:*(p::Polynomial, c::Union{Integer,Rational}) = c * p

function Base.:/(p::Polynomial, c::Union{Integer,Rational})
  iszero(c) && throw(DivideError())
  Polynomial(p.coeffs ./ c)
end

function Base.:^(p::Polynomial, n::Integer)
  n >= 0 || throw(DomainError(n, "The polynomial exponent must be nonnegative."))
  Base.power_by_squaring(p, n)
end

Base.literal_pow(::typeof(^), p::Polynomial, ::Val{n}) where {n} = p^n

function Base.show(io::IO, p::Polynomial)
  iszero(p) && return print(io, "0")
  first_term = true
  for k in degree(p):-1:0
    c = coefficient(p, k)
    iszero(c) && continue
    if first_term
      c < 0 && print(io, "-")
    else
      print(io, c < 0 ? " - " : " + ")
    end
    if k == 0 || abs(c) != 1
      magnitude = abs(c)
      print(io, denominator(magnitude) == 1 ? numerator(magnitude) : magnitude)
      k > 0 && print(io, "*")
    end
    k > 0 && print(io, "t")
    k > 1 && print(io, "^", k)
    first_term = false
  end
  nothing
end
