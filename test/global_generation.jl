@testset "Global generation" begin
  for X in (projective_space(3), Gr(2, 5), OGr(2, 7), SGr(2, 6),
    cayley_plane(), partial_flag_variety(TypeG2, 1))
    @test is_globally_generated(structure_sheaf(X))
    @test is_globally_generated(zero_bundle(X))
    @test is_globally_generated(anticanonical_bundle(X))
    @test !is_globally_generated(canonical_bundle(X))
  end

  X = Gr(2, 5)
  U, Q = universal_subbundle(X), universal_quotient_bundle(X)
  @test !is_globally_generated(U)
  @test is_globally_generated(dual(U))
  @test is_globally_generated(Q)
  @test is_globally_generated(symmetric_power(dual(U), 2))
  @test is_globally_generated(direct_sum(dual(U), Q))
  @test !is_globally_generated(direct_sum(dual(U), U))
  @test !is_globally_generated(direct_sum(zero_bundle(X), U))

  F = flag_variety(4, [1, 3])
  @test is_globally_generated(line_bundle(F, [0, 1]))
  @test !is_globally_generated(line_bundle(F, [1, -1]))
  @test is_globally_generated(line_bundle(F, [0, 0]))

  Y = projective_space(2)
  @test is_globally_generated(external_tensor_product(dual(U), line_bundle(Y, 1)))
  @test !is_globally_generated(external_tensor_product(dual(U), line_bundle(Y, -1)))

  E = symmetric_power(dual(U), 2)
  @test zero_locus(E; check_global_generation=true) == zero_locus(E)
  @test_throws ArgumentError zero_locus(U; check_global_generation=true)
  @test defining_bundle(zero_locus(U)) == U
  @test defining_bundle(zero_locus(U; check_global_generation=false)) == U
end
