# ═══════════════════════════════════════════════════════════════════════════════
#  External products of ambient bundles
# ═══════════════════════════════════════════════════════════════════════════════

export external_direct_sum, external_tensor_product

const _AmbientBundle = Union{CompletelyReducibleBundle,FilteredBundle}

# Lift every graded piece of `F` into one factor of `product_ambient`.
function _lift_bundle_to_product(
  product_ambient::PartialFlagVariety, F::FilteredBundle, offset::Int
)
  FilteredBundle(
    product_ambient,
    [
      _lift_bundle_to_product(product_ambient, piece, offset) for
      piece in graded_pieces(F)
    ],
  )
end

# Lift two ambient bundles to the common product ambient, placing the right
# bundle after the left Dynkin block.
function _lift_external_factors(
  product_ambient::PartialFlagVariety,
  E::_AmbientBundle,
  F::_AmbientBundle,
  right_offset::Int,
)
  (
    _lift_bundle_to_product(product_ambient, E, 0),
    _lift_bundle_to_product(product_ambient, F, right_offset),
  )
end

"""
    external_tensor_product(E, F) -> Bundle

Return the external tensor product
``p_X^*\\mathcal{E} \\otimes p_Y^*\\mathcal{F}`` on ``X \\times Y``, where
`E` lives on `X` and `F` lives on `Y`.

The arguments may be [`CompletelyReducibleBundle`](@ref)s or
[`FilteredBundle`](@ref)s. If either argument is filtered, the result carries
the induced tensor-product filtration.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X, Y = projective_space(1), projective_space(2);

julia> L = external_tensor_product(line_bundle(X, 2), line_bundle(Y, 3));

julia> L == line_bundle(product(X, Y), [2, 3])
true
```
"""
function external_tensor_product(E::_AmbientBundle, F::_AmbientBundle)
  product_ambient = product(variety(E), variety(F))
  right_offset = rank(dynkin_type(variety(E)))
  lifted_E, lifted_F = _lift_external_factors(product_ambient, E, F, right_offset)
  tensor_product(lifted_E, lifted_F)
end

"""
    external_direct_sum(
      E::CompletelyReducibleBundle,
      F::CompletelyReducibleBundle,
    ) -> CompletelyReducibleBundle

Return the external direct sum
``p_X^*\\mathcal{E} \\oplus p_Y^*\\mathcal{F}`` on ``X \\times Y``, where
`E` lives on `X` and `F` lives on `Y`.

This method is restricted to completely reducible bundles. A
[`FilteredBundle`](@ref) records only an ordered list of graded pieces, without
filtration indices that could canonically align the two summand filtrations.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X, Y = projective_space(1), projective_space(2);

julia> E = external_direct_sum(line_bundle(X, 2), line_bundle(Y, 3));

julia> E == direct_sum(line_bundle(product(X, Y), [2, 0]),
                       line_bundle(product(X, Y), [0, 3]))
true
```
"""
function external_direct_sum(
  E::CompletelyReducibleBundle, F::CompletelyReducibleBundle
)
  product_ambient = product(variety(E), variety(F))
  right_offset = rank(dynkin_type(variety(E)))
  lifted_E, lifted_F = _lift_external_factors(product_ambient, E, F, right_offset)
  direct_sum(lifted_E, lifted_F)
end
