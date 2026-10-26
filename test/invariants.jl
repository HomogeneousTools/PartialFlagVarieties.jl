@testset "Anticanonical volume and invariant reports" begin
  for n in (1, 2, 3, 21)
    @test volume(projective_space(n)) == big(n + 1)^n
  end
  @test volume(quadric(3)) == 54
  @test volume(full_flag_variety(TypeA{2})) == 48
  @test volume(projective_space(1) * projective_space(2)) == 54
  @test anticanonical_volume(projective_space(2)) == 9
  @test volume(projective_space(2)) isa BigInt

  X = projective_space(4)
  cubic = zero_locus(line_bundle(X, 3))
  @test volume(cubic) == 24
  @test volume(zero_locus(direct_sum(line_bundle(X, 2), line_bundle(X, 2)))) == 4
  @test volume(zero_locus(structure_sheaf(X))) == 0
  @test_throws ArgumentError volume(zero_locus(line_bundle(X, 5)))
  @test_throws ArgumentError volume(zero_locus(line_bundle(X, 6)))

  point = partial_flag_variety(TypeA{1}, Int[])
  @test volume(point) == 1
  Ipoint = invariants(point)
  @test Ipoint.dimension == 0
  @test isempty(Ipoint.chi_tangent)
  @test Ipoint.chi_cotangent == [1]

  I = invariants(projective_space(2); max_tangent_power=2)
  @test I.dimension == 2
  @test I.volume == 9
  @test I.chi_structure == 1
  @test I.chi_tangent == [8, 10]
  @test I.chi_cotangent == [1, -1]
  @test I.hodge_numbers == BigInt[1 0 0; 0 1 0; 0 0 1]
  @test invariants(projective_space(2); hodge=false).hodge_numbers === nothing
  @test isempty(invariants(projective_space(2); max_tangent_power=0).chi_tangent)
  @test_throws ArgumentError invariants(projective_space(2); max_tangent_power=-1)
  @test_throws ArgumentError invariants(projective_space(2); max_tangent_power=3)

  K3 = zero_locus(line_bundle(projective_space(3), 4))
  I = invariants(K3; hodge=false, max_tangent_power=2)
  @test I.volume === nothing
  @test I.chi_structure == 2
  @test I.chi_tangent == [-20, 2]
  @test I.chi_cotangent == [2, -20]
  @test I.hodge_numbers === nothing

  I = invariants(cubic; hodge=false)
  @test I.volume == 24
  @test I.chi_structure == 1
  @test I.chi_tangent == [-10]
  @test I.chi_cotangent == [1, 4]
end
