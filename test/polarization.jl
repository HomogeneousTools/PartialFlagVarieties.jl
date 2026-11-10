@testset "Stored polarizations" begin
  X = projective_space(2)
  F = structure_sheaf(X)
  original_hash = hash(X)
  @test polarization(X) == line_bundle(X, 1)
  @test set_polarization!(X, 2) === X
  @test hilbert_polynomial(F) == hilbert_polynomial(F, line_bundle(X, 2))
  @test hilbert_polynomial(F)(3) == euler_characteristic(line_bundle(X, 6))
  @test hilbert_polynomial(F, line_bundle(X, 1))(3) == 10
  @test polarization(X) == line_bundle(X, 2)
  @test X == projective_space(2)
  @test hash(X) == original_hash
  @test set_polarization!(X, nothing) === X
  @test polarization(X) == line_bundle(X, 1)

  Y = partial_flag_variety("A2", [1, 2]; polarization=[2, 3])
  @test picard_degrees(polarization(Y)) == [2, 3]
  @test picard_degrees(polarization(X * Y)) == [1, 2, 3]
  @test_throws ArgumentError set_polarization!(Y, 1)
  @test_throws ArgumentError set_polarization!(Y, [1, 0])
  @test_throws ArgumentError set_polarization!(X, structure_sheaf(X))
  @test_throws ArgumentError set_polarization!(X, tangent_bundle(X))
  @test_throws ArgumentError set_polarization!(X, line_bundle(projective_space(1), 1))

  Z = zero_locus(line_bundle(X, 2))
  @test polarization(Z) == line_bundle(Z, 1)
  set_polarization!(X, 3)
  @test polarization(Z) == line_bundle(Z, 3)
  @test hilbert_polynomial(Z)(2) == 13
  @test set_polarization!(Z, line_bundle(Z, 2)) === Z
  set_polarization!(X, 4)
  @test polarization(Z) == line_bundle(Z, 2)
  @test hilbert_polynomial(Z)(2) == 9
  @test hilbert_polynomial(tangent_bundle(Z))(2) == 11
  @test hilbert_polynomial(Z, polarization(Z)) == hilbert_polynomial(Z)
  @test hilbert_polynomial(tangent_bundle(Z), polarization(Z)) ==
    hilbert_polynomial(tangent_bundle(Z))
  @test hilbert_polynomial(Z, line_bundle(X, 1))(2) == 5
  @test Z == zero_locus(line_bundle(projective_space(2), 2))
  set_polarization!(Z, nothing)
  @test polarization(Z) == line_bundle(Z, 4)

  W = zero_locus(line_bundle(X, 2); polarization=2)
  @test polarization(W) == line_bundle(W, 2)
  @test polarization(Z * W) == line_bundle(Z * W, [4, 2])
  @test_throws ArgumentError set_polarization!(W, line_bundle(X, -1))
  @test_throws ArgumentError set_polarization!(
    W, line_bundle(zero_locus(line_bundle(X, 3)), 1)
  )
  @test_throws ArgumentError set_polarization!(W, cotangent_bundle(W))

  flag = full_flag_variety(TypeA{2}; polarization=[2, 3])
  @test hilbert_polynomial(tangent_bundle(flag)) ==
    hilbert_polynomial(tangent_bundle(flag), line_bundle(flag, [2, 3]))
end
