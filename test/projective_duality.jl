@testset "Projective duality" begin
  function check_dual(X, degree, defect=0, L=nothing)
    args = L === nothing ? (X,) : (X, L)
    @test projective_dual_degree(args...) == BigInt(degree)
    @test projective_dual_codimension(args...) == defect + 1
    @test projective_dual_defect(args...) == defect
  end

  @testset "Tevelev examples (arXiv numbering)" begin
    # Example 7.9: Veronese embeddings and classical discriminants.
    for n in 1:4, d in 2:4
      X = projective_space(n)
      check_dual(X, (n + 1) * (d - 1)^n, 0, line_bundle(X, d))
    end

    # Examples 7.6--7.8: plane curves with the degree-r polarization,
    # hypersurfaces (including surfaces), and a genus-one complete intersection.
    P2 = projective_space(2)
    for d in 2:5, r in 1:3
      C = zero_locus(line_bundle(P2, d))
      check_dual(C, d * (d + 2r - 3), 0, line_bundle(C, r))
    end
    for n in 1:3, d in 2:4
      P = projective_space(n + 1)
      check_dual(zero_locus(line_bundle(P, d)), d * (d - 1)^n)
    end
    P3 = projective_space(3)
    elliptic = zero_locus(2 * line_bundle(P3, 2))
    check_dual(elliptic, 8)
    check_dual(elliptic, 16, 0, line_bundle(P3, 2))
    K3 = zero_locus(line_bundle(P3, 4))
    check_dual(K3, 72, 0, line_bundle(K3, 2))

    # Example 2.18: determinantal duals of Segre varieties, including
    # rectangular formats with positive defect. See also Theorem 6.8.
    for m in 1:3, n in m:4
      check_dual(projective_space(m) * projective_space(n), binomial(n + 1, m), n - m)
    end
    check_dual(projective_space(1) * projective_space(1) * projective_space(1), 4)
    check_dual(projective_space(1) * projective_space(1) * projective_space(2), 6)

    # Examples 2.19--2.22: quadrics, Pfaffians, the self-dual spinor tenfold,
    # and the Cartan cubic dual to the Cayley plane.
    # Q^1 is represented as P^1: its minimal polarization is O(1),
    # so specify O(2) to get the conic embedding.
    Q1 = quadric(1)
    check_dual(Q1, 2, 0, line_bundle(Q1, 2))
    for n in 2:6
      check_dual(quadric(n), 2)
    end
    for n in (4, 6, 8)
      check_dual(Gr(2, n), n ÷ 2)
    end
    check_dual(Gr(2, 5), 5, 2)
    check_dual(Gr(2, 7), 14, 2)
    check_dual(OGr(5, 10), 12, 4)
    check_dual(cayley_plane(), 3)

    # Example 7.18: full flags of SL(3), with two independent Picard degrees.
    flag = full_flag_variety(TypeA{2})
    for (a, b) in ((1, 1), (1, 2), (2, 3), (3, 3))
      degree = 12 * (a^2 * b + a * b^2) - 6 * (a^2 + 4a*b + b^2) + 12 * (a+b) - 6
      check_dual(flag, degree, 0, line_bundle(flag, [a, b]))
    end
  end

  @testset "Issue #38 examples" begin
    for (X, degree) in ((Gr(3, 9), 120), (Gr(4, 8), 126),
      (OGr(8, 16), 240), (LGr(4, 8), 72))
      check_dual(X, degree)
    end
    X = projective_space(4)
    check_dual(zero_locus(line_bundle(X, 5)), 320)
    check_dual(Gr(3, 6), 4)
    check_dual(LGr(3, 6), 4)

    # General sections lower positive defect one step at a time (Theorem 1.30).
    X = Gr(2, 5)
    L = line_bundle(X, 1)
    check_dual(zero_locus(L), 5, 1)
    check_dual(zero_locus(L + L), 5, 0)
  end

  @testset "Defaults, cache reuse, and validation" begin
    PartialFlagVarieties.clear_caches!()
    X = projective_space(2)
    @test projective_dual_degree(X) === nothing
    check_dual(X, 3, 0, line_bundle(X, 2))
    after_degree = PartialFlagVarieties.cache_info().projective_dual
    @test projective_dual_codimension(X, line_bundle(X, 2)) == 1
    @test PartialFlagVarieties.cache_info().projective_dual.hits == after_degree.hits + 1
    @test projective_dual_degree(X) === nothing
    set_polarization!(X, 3)
    check_dual(X, 12)
    @test projective_dual_degree(X, line_bundle(X, 2)) == 3
    @test polarization(X) == line_bundle(X, 3)
    @test projective_dual_degree(projective_space(2)) === nothing
    set_polarization!(X, nothing)
    @test projective_dual_degree(X) === nothing

    Z = zero_locus(line_bundle(X, 2))
    check_dual(Z, 2)
    set_polarization!(X, 2)
    check_dual(Z, 6)
    set_polarization!(Z, 3)
    check_dual(Z, 10)
    check_dual(Z, 2, 0, line_bundle(Z, 1))
    @test projective_dual_degree(Z, polarization(Z)) == 10
    set_polarization!(Z, nothing)
    check_dual(Z, 6)

    @test_throws ArgumentError projective_dual_degree(X, tangent_bundle(X))
    @test_throws ArgumentError projective_dual_degree(
      X, line_bundle(projective_space(3), 1)
    )
    @test_throws ArgumentError projective_dual_degree(X, structure_sheaf(X))
    @test_throws ArgumentError projective_dual_degree(X, line_bundle(X, -1))
    Y = projective_space(1) * projective_space(1)
    @test_throws ArgumentError projective_dual_degree(Y, line_bundle(Y, [1, 0]))
    @test_throws ArgumentError projective_dual_degree(
      Z, line_bundle(zero_locus(line_bundle(X, 3)), 1)
    )
    @test_throws ArgumentError projective_dual_degree(Z, cotangent_bundle(Z))
    PartialFlagVarieties.configure_caches!(; budget=64 * 1024^2)
    @test PartialFlagVarieties.cache_info().projective_dual.maxsize > 0
    PartialFlagVarieties.clear_caches!()
    @test PartialFlagVarieties.cache_info().projective_dual.currentsize == 0
  end

  @testset "Linear spaces, empty loci, and points" begin
    P3 = projective_space(3)
    sources = (projective_space(1), P3, zero_locus(line_bundle(P3, 1)),
      zero_locus(3 * line_bundle(P3, 1)), zero_locus(structure_sheaf(P3)))
    for X in sources
      @test projective_dual_degree(X) === nothing
      @test projective_dual_codimension(X) === nothing
      @test projective_dual_defect(X) === nothing
    end
    P1 = projective_space(1)
    for m in 2:4
      check_dual(zero_locus(line_bundle(P1, m)), m)
    end
  end
end
