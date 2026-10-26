@testset "Bounds for symbolic cohomology" begin
  PFV = PartialFlagVarieties
  X = SGr(2, 6)
  F = tensor_product(exterior_power(filtered_cotangent_bundle(X), 3), line_bundle(X, -1))
  result = cohomology(F; with_bounds=true)
  @test result.values.entries == cohomology(F).entries
  @test result.values[4] == AffineExpr(2) - symbolic_variable(0)
  @test result.bounds[0] == (lower=0, upper=1)

  twisted = twisted_hodge_numbers(X, -1; with_bounds=true)
  @test twisted.values == twisted_hodge_numbers(X, -1)
  @test twisted.bounds[0] == (lower=0, upper=1)

  flag = full_flag_variety(TypeB{2})
  polyvectors = hochschild_cohomology(flag; with_bounds=true)
  @test polyvectors.values.data == hochschild_cohomology(flag).data
  @test !isempty(polyvectors.bounds)
  ids = Set(id for entry in polyvectors.values.data for id in keys(entry.coeffs))
  @test Set(keys(polyvectors.bounds)) == ids

  # The upper bound comes from a Koszul map rank, not just nonnegativity
  # of the two final cohomology entries.
  P4 = projective_space(4)
  cubic = zero_locus(line_bundle(P4, 3))
  result = cohomology(tangent_bundle(cubic); with_bounds=true)
  @test result.values.entries == cohomology(tangent_bundle(cubic)).entries
  @test result.bounds[0] == (lower=10, upper=34)

  X = Gr(2, 7)
  E = symmetric_power(dual(universal_subbundle(X)), 2)
  Z = zero_locus(direct_sum(E, E))
  result = hodge_numbers(Z; with_bounds=true)
  @test result.values == hodge_numbers(Z)
  @test !isempty(result.bounds)
  ids = Set(id for entry in result.values for id in keys(entry.coeffs))
  @test Set(keys(result.bounds)) == ids
  @test all(
    interval ->
      interval.lower >= 0 &&
      (interval.upper === nothing || interval.lower <= interval.upper),
    values(result.bounds))

  P2 = projective_space(2)
  conic = zero_locus(line_bundle(P2, 2))
  @test isempty(cohomology(structure_sheaf(P2); with_bounds=true).bounds)
  @test isempty(cohomology(structure_sheaf(P2); characters=true, with_bounds=true).bounds)
  @test isempty(hodge_numbers(P2; with_bounds=true).bounds)
  @test isempty(hodge_numbers(conic; with_bounds=true).bounds)
  @test isempty(hodge_numbers_symbolic(conic; with_bounds=true).bounds)
  @test isempty(hochschild_cohomology(conic; with_bounds=true).bounds)
  @test isempty(hochschild_cohomology(P2; with_bounds=true).bounds)
  @test isempty(twisted_hodge_numbers(P2, 0; with_bounds=true).bounds)
  @test isempty(twisted_hodge_numbers(P2, [0]; with_bounds=true).bounds)
  @test isempty(twisted_hodge_numbers(conic, 0; with_bounds=true).bounds)
  @test isempty(twisted_hodge_numbers(conic, [0]; with_bounds=true).bounds)

  @testset "Elimination, renumbering and unbounded parameters" begin
    result = PFV._with_cohomology_bounds(true) do
      x, y = symbolic_variable(7), symbolic_variable(20)
      PFV._record_nonnegative!(AffineExpr(8) - x)
      entries = AffineExpr[x]
      PFV._substitute_var!(entries, 7, AffineExpr(5) + y)
      PFV._renumber_variables!(entries)
    end
    @test result.values == [AffineExpr(5) + symbolic_variable(0)]
    @test result.bounds == Dict(0 => (lower=0, upper=3))

    result = PFV._with_cohomology_bounds(true) do
      PFV._renumber_variables!([symbolic_variable(9)])
    end
    @test result.bounds[0] == (lower=0, upper=nothing)

    x, y = symbolic_variable(0), symbolic_variable(1)
    intervals = PFV._variable_intervals([x - AffineExpr(2), AffineExpr(9) - x, x - 2y])
    @test intervals[0] == (lower=2, upper=9)
    @test intervals[1] == (lower=0, upper=4)
    @test_throws ErrorException PFV._variable_intervals([
      x - AffineExpr(3), AffineExpr(2) - x
    ])
  end

  @testset "Nested calls and task isolation" begin
    @test PFV._bounds_context() === nothing
    PFV._with_cohomology_bounds(true) do
      context = PFV._bounds_context()
      cohomology(F; with_bounds=true)
      @test PFV._bounds_context() === context
      @test_throws ArgumentError twisted_hodge_numbers(
        P2, tangent_bundle(P2); with_bounds=true
      )
      @test PFV._bounds_context() === context
      @test isempty(context.constraints)
    end
    @test PFV._bounds_context() === nothing
    tasks = [@async(cohomology(F; with_bounds=true)) for _ in 1:2]
    @test all(task -> fetch(task).bounds == Dict(0 => (lower=0, upper=1)), tasks)
    @test PFV._bounds_context() === nothing
  end
end
