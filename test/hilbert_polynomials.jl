@testset "Symbolic Weyl Hilbert polynomials" begin
  P21 = projective_space(21)
  polynomial = hilbert_polynomial(structure_sheaf(P21))
  @test length(polynomial) == 22
  @test last(polynomial) == 1 // factorial(big(21))
  @test evalpoly(big(21), polynomial) == binomial(big(42), 21)
  @test hilbert_polynomial(structure_sheaf(P21); max_degree=1) == polynomial

  # Compare with BWB, including negative twists, multiplicities, products,
  # exceptional types and both root lengths in B, C, F and G.
  varieties = (Gr(2, 5), OGr(2, 7), SGr(2, 6),
    partial_flag_variety(TypeF4, 4), partial_flag_variety(TypeG2, 1),
    cayley_plane(), projective_space(1) * projective_space(2))
  for X in varieties
    L = line_bundle(X, ones(Int, picard_rank(X)))
    F = direct_sum(cotangent_bundle(X), 2 * structure_sheaf(X))
    polynomial = hilbert_polynomial(F, L)
    weight = p_dominant_weight(only(components(L)))
    for t in (-3, -1, 0, 2, 4)
      Lt = CompletelyReducibleBundle(X, t * weight)
      @test evalpoly(big(t), polynomial) == euler_characteristic(tensor_product(F, Lt))
    end
    @test hilbert_polynomial(zero_bundle(X), L) == Rational{BigInt}[0]
  end

  X = SGr(2, 6)
  F = filtered_tangent_bundle(X)
  @test hilbert_polynomial(F) == hilbert_polynomial(total_bundle(F))
  @test hilbert_polynomial(FilteredBundle(X, CompletelyReducibleBundle[])) == [0]

  X = projective_space(2)
  Z = zero_locus(line_bundle(X, 2))
  @test hilbert_polynomial(Z) == Rational{BigInt}[1, 2]
  @test hilbert_polynomial(restrict(Z, line_bundle(X, 1))) == Rational{BigInt}[3, 2]
  @test hilbert_polynomial(tangent_bundle(Z)) == Rational{BigInt}[3, 2]
  @test hilbert_polynomial(zero_bundle(Z)) == Rational{BigInt}[0]
  @test hilbert_polynomial(zero_locus(structure_sheaf(X))) == Rational{BigInt}[0]
  @test hilbert_polynomial(zero_locus(zero_bundle(X))) ==
    hilbert_polynomial(structure_sheaf(X))

  P3 = projective_space(3)
  K3 = zero_locus(line_bundle(P3, 4))
  @test hilbert_polynomial(K3) == Rational{BigInt}[2, 0, 2]
  @test hilbert_polynomial(tangent_bundle(K3)) == Rational{BigInt}[-20, 0, 4]

  X = projective_space(4)
  Z = zero_locus(line_bundle(X, 5))
  L = line_bundle(X, 1)
  F = tensor_product(tangent_bundle(Z), cotangent_bundle(Z))
  polynomial = hilbert_polynomial(F, L)
  for t in (-2, 0, 1, 3)
    @test evalpoly(big(t), polynomial) ==
      euler_characteristic(tensor_product(F, line_bundle(Z, t)))
  end
  @test_throws ArgumentError hilbert_polynomial(F, tangent_bundle(X))
  @test_throws ArgumentError hilbert_polynomial(F, line_bundle(P3, 1))
  @test_throws ArgumentError hilbert_polynomial(Z, line_bundle(P3, 1))
  @test_throws ArgumentError hilbert_polynomial(
    filtered_tangent_bundle(X), tangent_bundle(X)
  )
end
