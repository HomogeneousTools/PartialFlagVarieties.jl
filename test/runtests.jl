using Test
using PartialFlagVarieties
using Lie
using StaticArrays

@testset "PartialFlagVarieties.jl" begin

  # ═══════════════════════════════════════════════════════════════════════════
  #  Lie.jl extensions: cartan_type, parse_dynkin_type
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "cartan_type" begin
    # A3 Cartan matrix
    C_A3 = cartan_matrix(TypeA{3})
    @test cartan_type(C_A3) == [(:A, 3)]
    ct_ord = cartan_type_with_ordering(C_A3)
    @test ct_ord[1] == [(:A, 3)]
    @test length(ct_ord[2]) == 3

    # B3 Cartan matrix
    C_B3 = cartan_matrix(TypeB{3})
    @test cartan_type(C_B3) == [(:B, 3)]

    # Product type A1 × A2 (block-diagonal Cartan matrix)
    C_prod = cartan_matrix(ProductDynkinType{Tuple{TypeA{1},TypeA{2}}})
    ct = cartan_type(C_prod)
    @test length(ct) == 2
    # Should contain A1 and A2 components (order may vary)
    types = Set(ct)
    @test (:A, 1) in types
    @test (:A, 2) in types
  end

  @testset "parse_dynkin_type" begin
    @test parse_dynkin_type("A3") === TypeA{3}
    @test parse_dynkin_type("B4") === TypeB{4}
    @test parse_dynkin_type("G2") === TypeG2
    @test parse_dynkin_type("E6") === TypeE{6}

    # Product types
    DT = parse_dynkin_type("A2xB3")
    @test DT === ProductDynkinType{Tuple{TypeA{2}, TypeB{3}}}

    # Whitespace tolerance
    @test parse_dynkin_type(" A3 ") === TypeA{3}
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  MarkedDynkinType
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "MarkedDynkinType basics" begin
    MDT = MarkedDynkinType{TypeA{4},(2,)}

    @test marked_nodes(MDT) == (2,)
    @test unmarked_nodes(MDT) == (1, 3, 4)
    @test central_rank(MDT) == 1
    @test levi_rank(MDT) == 3
    @test rank(MDT) == 4

    # Full flag
    MDT_full = MarkedDynkinType{TypeA{3}, (1, 2, 3)}
    @test is_full_flag(MDT_full)
    @test levi_type(MDT_full) === nothing
    @test central_rank(MDT_full) == 3
    @test levi_rank(MDT_full) == 0

    # Two marked nodes
    MDT2 = MarkedDynkinType{TypeA{4}, (1, 3)}
    @test marked_nodes(MDT2) == (1, 3)
    @test unmarked_nodes(MDT2) == (2, 4)
    @test central_rank(MDT2) == 2
  end

  @testset "MarkedDynkinType constructors" begin
    mdt1 = MarkedDynkinType(TypeA{3}, (2,))
    @test mdt1 isa MarkedDynkinType{TypeA{3},(2,)}

    mdt2 = MarkedDynkinType(TypeB{4}, [1, 3])
    @test mdt2 isa MarkedDynkinType{TypeB{4},(1,3)}

    mdt3 = MarkedDynkinType(TypeD{5}, 5)
    @test mdt3 isa MarkedDynkinType{TypeD{5},(5,)}
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Levi type identification
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Levi type computation" begin
    # Gr(2,5) = A4/P2 → Levi = A1 × A2
    @test levi_type(MarkedDynkinType{TypeA{4},(2,)}) ==
          ProductDynkinType{Tuple{TypeA{1},TypeA{2}}}

    # Gr(1,5) = A4/P1 → Levi = A3
    @test levi_type(MarkedDynkinType{TypeA{4},(1,)}) == TypeA{3}

    # D5/P5 → Levi = A4
    @test levi_type(MarkedDynkinType{TypeD{5},(5,)}) == TypeA{4}

    # B3/P1 → Levi = B2
    @test levi_type(MarkedDynkinType{TypeB{3},(1,)}) == TypeB{2}

    # E6/P1 → Levi = D5
    @test levi_type(MarkedDynkinType{TypeE{6},(1,)}) == TypeD{5}

    # E6/P2 → Levi = A4 × ... (compute and check rank)
    @test levi_rank(MarkedDynkinType{TypeE{6},(2,)}) == 5
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  PartialFlagVariety constructors
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "PartialFlagVariety constructors" begin
    V = partial_flag_variety(TypeA{3}, (2,))
    @test V isa PartialFlagVariety
    @test marked_nodes(V) == (2,)
    @test dynkin_type(V) == TypeA{3}

    V_full = full_flag_variety(TypeA{2})
    @test is_full_flag(V_full)

    # String constructor
    V_str = PartialFlagVariety("A3", [2])
    @test dimension(V_str) == 4
    @test marked_nodes(V_str) == (2,)

    # Product type string constructor
    V_prod = PartialFlagVariety("A2xB3", [1, 4])
    @test V_prod isa PartialFlagVariety
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Dimensions (on PartialFlagVariety)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Dimensions" begin
    # Projective spaces: dim(ℙⁿ) = n
    @test dimension(partial_flag_variety(TypeA{1}, (1,))) == 1
    @test dimension(partial_flag_variety(TypeA{2}, (1,))) == 2
    @test dimension(partial_flag_variety(TypeA{4}, (1,))) == 4

    # Grassmannians: dim(Gr(k,n)) = k(n-k)
    @test dimension(Gr(2, 4)) == 4
    @test dimension(Gr(2, 5)) == 6
    @test dimension(Gr(3, 6)) == 9

    # Full flag: dim(G/B) = |Φ⁺|
    @test dimension(partial_flag_variety(TypeA{2}, (1, 2))) == 3
    @test dimension(full_flag_variety(TypeA{3})) == 6

    # Exceptional: Cayley plane OP² = E6/P1 has dim 16
    @test dimension(cayley_plane()) == 16

    # Quadrics: Q_n has dim n
    @test dimension(quadric(3)) == 3
    @test dimension(quadric(4)) == 4

    # OGr(5,10) = D5/P5: spinor variety, dim 10
    @test dimension(OGr(5, 10)) == 10
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Picard rank
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Picard rank" begin
    @test picard_rank(Gr(2, 5)) == 1
    @test picard_rank(partial_flag_variety(TypeA{3}, (1, 3))) == 2
    @test picard_rank(full_flag_variety(TypeA{3})) == 3
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Euler characteristics
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Euler characteristics" begin
    # χ(ℙⁿ) = n + 1
    @test euler_characteristic(projective_space(4)) == 5

    # χ(Gr(k,n)) = C(n,k)
    @test euler_characteristic(Gr(2, 4)) == 6
    @test euler_characteristic(Gr(2, 5)) == 10
    @test euler_characteristic(Gr(2, 6)) == 15
    @test euler_characteristic(Gr(3, 6)) == 20

    # χ(OP²) = 27 (Cayley plane)
    @test euler_characteristic(cayley_plane()) == 27

    # χ(E7/P7) = 56 (Freudenthal variety)
    @test euler_characteristic(freudenthal_variety()) == 56

    # χ(OGr(5,10)) = 2⁴ = 16
    @test euler_characteristic(OGr(5, 10)) == 16

    # χ(G/B) = |W|
    @test euler_characteristic(full_flag_variety(TypeA{2})) == 6
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Betti numbers
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Betti numbers" begin
    # ℙⁿ: all Betti numbers = 1
    @test betti_numbers(projective_space(2)) == [1, 1, 1]
    @test betti_numbers(projective_space(4)) == [1, 1, 1, 1, 1]

    # Gr(2,4): Betti = 1, 1, 2, 1, 1
    @test betti_numbers(Gr(2, 4)) == [1, 1, 2, 1, 1]

    # Gr(2,5): Betti = 1, 1, 2, 2, 2, 1, 1
    @test betti_numbers(Gr(2, 5)) == [1, 1, 2, 2, 2, 1, 1]

    # Sum of Betti = χ
    @test sum(betti_numbers(cayley_plane())) == 27
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Decomposition matrix
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Decomposition matrix" begin
    MDT = MarkedDynkinType{TypeA{3},(2,)}
    M = decomposition_matrix(MDT)
    Minv = decomposition_matrix_inv(MDT)

    # M * Minv = I
    @test M * Minv ≈ SMatrix{3,3,Rational{Int}}(1, 0, 0, 0, 1, 0, 0, 0, 1)

    # Unmarked rows are identity
    @test M[1, 1] == 1 && M[1, 2] == 0 && M[1, 3] == 0
    @test M[3, 3] == 1
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Classification predicates
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Classification predicates" begin
    @test is_generalized_grassmannian(Gr(2, 4))
    @test !is_generalized_grassmannian(partial_flag_variety(TypeA{3}, (1, 3)))

    @test is_full_flag(full_flag_variety(TypeA{2}))
    @test !is_full_flag(Gr(2, 4))

    # Cominuscule
    @test is_cominuscule(Gr(2, 5))
    @test is_cominuscule(quadric(8))                            # D5/P1
    @test is_cominuscule(OGr(5, 10))                            # D5/P5
    @test is_cominuscule(cayley_plane())                        # E6/P1
    @test !is_cominuscule(partial_flag_variety(TypeB{3}, (2,)))

    # Minuscule
    @test is_minuscule(Gr(2, 5))
    @test is_minuscule(partial_flag_variety(TypeB{3}, (3,)))
    @test is_minuscule(OGr(5, 10))
    @test !is_minuscule(partial_flag_variety(TypeB{3}, (1,)))

    # Adjoint
    @test is_adjoint(partial_flag_variety(TypeA{3}, (1, 3)))    # ℙ(T*ℙ³)
    @test is_adjoint(adjoint_variety(TypeB{3}))
    @test is_adjoint(adjoint_variety(TypeG2))
    @test is_adjoint(partial_flag_variety(TypeE{6}, (2,)))

    # Coadjoint
    @test is_coadjoint(coadjoint_variety(TypeB{3}))
    @test is_coadjoint(coadjoint_variety(TypeG2))
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Tangent weights
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Tangent weights" begin
    # ℙⁿ: tangent bundle is irreducible → 1 tangent weight
    @test length(tangent_weights(MarkedDynkinType{TypeA{4},(1,)})) == 1

    # Gr(2,4) = A3/P2: T = S* ⊗ Q (irreducible under Levi) → 1 tangent weight
    @test length(tangent_weights(MarkedDynkinType{TypeA{3},(2,)})) == 1

    # Full flag: each positive root is maximal → #tangent weights = dim
    @test length(tangent_weights(MarkedDynkinType{TypeA{2},(1,2)})) == 3
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Roots decomposition
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Root decomposition" begin
    MDT = MarkedDynkinType{TypeA{3},(2,)}

    nonpar = positive_nonparabolic_roots(MDT)
    par = positive_parabolic_roots(MDT)

    # Total = n_positive_roots(A3) = 6
    @test length(nonpar) + length(par) == Lie.n_positive_roots(TypeA{3})

    # dim(G/P) = number of nonparabolic positive roots
    X = partial_flag_variety(TypeA{3}, (2,))
    @test length(nonpar) == dimension(X)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  IrrepLevi
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "IrrepLevi construction" begin
    MDT = MarkedDynkinType{TypeA{3},(2,)}

    # Fundamental weight ω₁
    ω₁ = fundamental_weight(TypeA{3}, 1)
    rep = IrrepLevi(MDT, ω₁)
    @test length(central_part(rep)) == 1
    @test semisimple_part(rep) isa WeightLatticeElem

    # Round-trip: to_ambient → IrrepLevi → to_ambient
    λ_back = to_ambient_weight(MDT, rep)
    @test λ_back == ω₁
  end

  @testset "IrrepLevi round-trip" begin
    MDT = MarkedDynkinType{TypeA{4},(2,)}

    for i in 1:4
      ω = fundamental_weight(TypeA{4}, i)
      rep = IrrepLevi(MDT, ω)
      @test to_ambient_weight(MDT, rep) == ω
    end
  end

  @testset "IrrepLevi dual" begin
    MDT = MarkedDynkinType{TypeA{3},(2,)}
    ω₁ = fundamental_weight(TypeA{3}, 1)
    rep = IrrepLevi(MDT, ω₁)
    d = dual(rep)
    # Central part is negated under dual
    @test central_part(d) == -central_part(rep)
    # P-dominant weight is stored and accessible
    @test p_dominant_weight(rep) == ω₁
    # Round-trip: dual of dual recovers original
    @test p_dominant_weight(dual(d)) == ω₁
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  CompletelyReducibleBundle
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Structure sheaf" begin
    X = Gr(2, 4)
    O = structure_sheaf(X)
    @test rank_bundle(O) == 1
    @test n_components(O) == 1
    @test variety(O) === X
  end

  @testset "Line bundles" begin
    X = projective_space(4)
    L = line_bundle(X, 1)
    @test rank_bundle(L) == 1
    @test variety(L) === X
  end

  @testset "Line bundle: Picard rank check" begin
    # Picard rank 1 → single-int form works
    X1 = Gr(2, 5)
    @test picard_rank(X1) == 1
    L = line_bundle(X1, 2)
    @test rank_bundle(L) == 1

    # Picard rank > 1 → single-int form errors
    X2 = partial_flag_variety(TypeA{3}, (1, 3))
    @test picard_rank(X2) == 2
    @test_throws ArgumentError line_bundle(X2, 1)

    # Multi-node form works for Picard rank > 1
    L2 = line_bundle(X2, [2, 1])
    @test rank_bundle(L2) == 1

    # Multi-node form also works for Picard rank 1
    L3 = line_bundle(X1, [3])
    @test rank_bundle(L3) == 1

    # Wrong number of degrees errors
    @test_throws ArgumentError line_bundle(X2, [1, 2, 3])
  end

  @testset "Tangent and cotangent bundles" begin
    # Gr(2,4): T = S* ⊗ Q, both have rank 2, so T has rank 4
    X = Gr(2, 4)
    T = tangent_bundle(X)
    Ω = cotangent_bundle(X)
    @test rank_bundle(T) == 4
    @test rank_bundle(Ω) == 4
    @test variety(T) === X

    # ℙ⁴: T has rank 4
    X2 = projective_space(4)
    T2 = tangent_bundle(X2)
    @test rank_bundle(T2) == 4
  end

  @testset "Bundle operations" begin
    X = projective_space(4)
    O = structure_sheaf(X)
    T = tangent_bundle(X)

    # O ⊗ T = T
    @test rank_bundle(tensor_product(O, T)) == rank_bundle(T)

    # ⊕
    @test rank_bundle(direct_sum(T, O)) == rank_bundle(T) + 1

    # ⋀⁰T = O, ⋀¹T = T
    @test rank_bundle(exterior_power(T, 0)) == 1
    @test rank_bundle(exterior_power(T, 1)) == rank_bundle(T)

    # variety is propagated through operations
    @test variety(tensor_product(O, T)) === X
    @test variety(direct_sum(T, O)) === X
    @test variety(exterior_power(T, 2)) === X
  end

  @testset "Twist" begin
    X = projective_space(4)
    O = structure_sheaf(X)
    Ot = twist(O, 1, 3)
    @test rank_bundle(Ot) == 1
    @test variety(Ot) === X
  end

  @testset "Determinant bundle" begin
    X = Gr(2, 4)
    T = tangent_bundle(X)
    det_T = det_bundle(T)
    @test rank_bundle(det_T) == 1
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Cohomology
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Cohomology: H*(ℙⁿ, 𝒪)" begin
    X = projective_space(4)
    O = structure_sheaf(X)
    H = cohomology(O)
    d = dimensions(H)

    # H⁰(ℙ⁴, 𝒪) = 1, Hⁱ = 0 for i > 0
    @test d[0] == 1
    for i in 1:4
      @test d[i] == 0
    end
  end

  @testset "Cohomology: H*(ℙⁿ, 𝒪(k))" begin
    X = projective_space(4)

    # H⁰(ℙ⁴, 𝒪(1)) = 5 (standard rep of SL(5))
    L1 = line_bundle(X, 1)
    H1 = dimensions(cohomology(L1))
    @test H1[0] == 5

    # H⁰(ℙ⁴, 𝒪(2)) = C(6,2) = 15
    L2 = twist(structure_sheaf(X), 1, 2)
    H2 = dimensions(cohomology(L2))
    @test H2[0] == 15

    # H⁰(ℙ⁴, 𝒪(3)) = C(7,3) = 35
    L3 = twist(structure_sheaf(X), 1, 3)
    H3 = dimensions(cohomology(L3))
    @test H3[0] == 35
  end

  @testset "Cohomology: Euler characteristic" begin
    X = projective_space(4)

    # χ(ℙ⁴, 𝒪) = 1
    @test euler_characteristic(structure_sheaf(X)) == 1

    # χ(ℜ⁴, 𝒪(1)) = 5
    @test euler_characteristic(line_bundle(X, 1)) == 5
  end

  @testset "Cohomology: 0-based indexing" begin
    X = projective_space(2)
    O = structure_sheaf(X)
    H = dimensions(cohomology(O))
    @test firstindex(H) == 0
    @test lastindex(H) == 2
    @test length(H) == 3
    @test H[0] == 1
    @test_throws BoundsError H[-1]
    @test_throws BoundsError H[3]
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Constructions
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Grassmannians" begin
    V = Gr(2, 5)
    @test dimension(V) == 6
    @test euler_characteristic(V) == 10

    V = Gr(3, 6)
    @test dimension(V) == 9
    @test euler_characteristic(V) == 20
  end

  @testset "Projective space" begin
    V = projective_space(4)
    @test dimension(V) == 4
    @test betti_numbers(V) == [1, 1, 1, 1, 1]
  end

  @testset "Flag varieties" begin
    V = flag_variety(4, (1, 2))
    @test dimension(V) == 5
  end

  @testset "Orthogonal Grassmannians" begin
    V = OGr(5, 10)
    @test dimension(V) == 10
    @test euler_characteristic(V) == 16

    V = OGr(1, 7)  # Q₅ = B₃/P₁
    @test dimension(V) == 5
  end

  @testset "Symplectic Grassmannians" begin
    V = SGr(2, 6)
    @test dimension(V) == 7

    V = LGr(3)
    @test dimension(V) == 6
  end

  @testset "Quadrics" begin
    V = quadric(4)
    @test dimension(V) == 4

    V = quadric(3)
    @test dimension(V) == 3
  end

  @testset "Exceptional varieties" begin
    V = cayley_plane()
    @test dimension(V) == 16
    @test euler_characteristic(V) == 27

    V = freudenthal_variety()
    @test dimension(V) == 27
    @test euler_characteristic(V) == 56
  end

  @testset "Adjoint varieties" begin
    V = adjoint_variety(TypeG2)
    @test dimension(V) == 5

    V = adjoint_variety(TypeB{3})
    @test dimension(V) == 7
  end

  @testset "Coadjoint varieties" begin
    V = coadjoint_variety(TypeG2)
    @test dimension(V) == 5

    V = coadjoint_variety(TypeB{3})
    @test dimension(V) == 5
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Marked Dynkin diagram display
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Marked Dynkin diagram" begin
    diagram = marked_dynkin_diagram(MarkedDynkinType{TypeA{4},(2,)})
    @test occursin("×", diagram)
    @test occursin("○", diagram)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Consistency checks
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Consistency: dim = #nonparabolic roots" begin
    for (DT, marks) in [
      (TypeA{4}, (2,)),
      (TypeB{3}, (1,)),
      (TypeD{5}, (5,)),
      (TypeE{6}, (1,)),
      (TypeG2, (1,)),
    ]
      MDT = MarkedDynkinType{DT,marks}
      X = partial_flag_variety(DT, marks)
      @test dimension(X) == length(positive_nonparabolic_roots(MDT))
    end
  end

  @testset "Consistency: χ = sum(betti)" begin
    for (DT, marks) in [
      (TypeA{3}, (2,)),
      (TypeA{4}, (1,)),
      (TypeB{3}, (1,)),
      (TypeD{5}, (5,)),
      (TypeE{6}, (1,)),
    ]
      X = partial_flag_variety(DT, marks)
      @test euler_characteristic(X) == sum(betti_numbers(X))
    end
  end

  @testset "Consistency: rank(T) = dim(X)" begin
    for (DT, marks) in [
      (TypeA{3}, (2,)),
      (TypeA{4}, (1,)),
      (TypeB{3}, (1,)),
    ]
      X = partial_flag_variety(DT, marks)
      @test Int(rank_bundle(tangent_bundle(X))) == dimension(X)
    end
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  levi_permutation
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "levi_permutation" begin
    # For type A, canonical ordering = natural ordering → identity permutation
    perm_A4 = levi_permutation(MarkedDynkinType{TypeA{4},(2,)})
    @test collect(perm_A4) == [1, 2, 3]

    # D4/P1: unmarked nodes (2,3,4) of D4 form an A3 sub-diagram
    # cartan_type_with_ordering finds the canonical A3 ordering = [2,1,3]
    perm_D4 = levi_permutation(MarkedDynkinType{TypeD{4},(1,)})
    @test collect(perm_D4) == [2, 1, 3]

    # B3/P1: unmarked nodes (2,3) form a B2 sub-diagram → identity
    perm_B3 = levi_permutation(MarkedDynkinType{TypeB{3},(1,)})
    @test length(perm_B3) == 2

    # G2/P1: unmarked node (2) → trivial A1 sub-diagram
    perm_G2 = levi_permutation(MarkedDynkinType{TypeG2,(1,)})
    @test collect(perm_G2) == [1]
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  IrrepLevi round-trip for non-A Dynkin types
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "IrrepLevi round-trip: D types" begin
    MDT = MarkedDynkinType{TypeD{4},(1,)}
    for i in 1:4
      ω = fundamental_weight(TypeD{4}, i)
      rep = IrrepLevi(MDT, ω)
      @test to_ambient_weight(MDT, rep) == ω
    end
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  IrrepLevi display: ambient node indices
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "IrrepLevi display uses ambient node indices" begin
    # Gr(3,8) = A7/P3: tangent bundle weight should show ω7 (not ω6)
    X = Gr(3, 8)
    T = tangent_bundle(X)
    s = sprint(show, T)
    @test occursin("ω7", s)
    @test !occursin("ω6", s)

    # D4/P1: tangent bundle should show ω2 (ambient D4 node 2)
    X_D4 = partial_flag_variety(TypeD{4}, (1,))
    T_D4 = tangent_bundle(X_D4)
    s_D4 = sprint(show, T_D4)
    @test occursin("ω2", s_D4)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Scalar multiplication on bundles (n * E = n-fold direct sum)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Scalar multiplication on bundles" begin
    X = Gr(2, 4)
    T = tangent_bundle(X)

    @test rank_bundle(2 * T) == 2 * rank_bundle(T)
    @test n_components(2 * T) == 2 * n_components(T)
    @test rank_bundle(T * 3) == 3 * rank_bundle(T)
    @test rank_bundle(1 * T) == rank_bundle(T)
    @test rank_bundle(0 * T) == 0
    @test variety(2 * T) === X
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Exterior power ranks (from reference output files)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Exterior power ranks: A-types" begin
    # A3/P2 = Gr(2,4), dim=4: ranks = binomial(4, k)
    T = tangent_bundle(Gr(2, 4))
    @test rank_bundle(exterior_power(T, 0)) == 1
    @test rank_bundle(exterior_power(T, 1)) == 4
    @test rank_bundle(exterior_power(T, 2)) == 6
    @test rank_bundle(exterior_power(T, 3)) == 4
    @test rank_bundle(exterior_power(T, 4)) == 1

    # A4/P2 = Gr(2,5), dim=6: ranks = binomial(6, k)
    T2 = tangent_bundle(Gr(2, 5))
    @test rank_bundle(exterior_power(T2, 0)) == 1
    @test rank_bundle(exterior_power(T2, 1)) == 6
    @test rank_bundle(exterior_power(T2, 2)) == 15
    @test rank_bundle(exterior_power(T2, 3)) == 20
    @test rank_bundle(exterior_power(T2, 4)) == 15
    @test rank_bundle(exterior_power(T2, 5)) == 6
    @test rank_bundle(exterior_power(T2, 6)) == 1

    # A7/P3 = Gr(3,8), dim=15: spot-check ranks = binomial(15, k)
    T3 = tangent_bundle(Gr(3, 8))
    @test rank_bundle(exterior_power(T3, 1)) == 15
    @test rank_bundle(exterior_power(T3, 2)) == 105
    @test rank_bundle(exterior_power(T3, 3)) == 455
  end

  @testset "Exterior power ranks: B/C/D/G types" begin
    # B3/P1, dim=5: ranks = binomial(5, k)
    T_B3 = tangent_bundle(partial_flag_variety(TypeB{3}, (1,)))
    @test rank_bundle(exterior_power(T_B3, 1)) == 5
    @test rank_bundle(exterior_power(T_B3, 2)) == 10
    @test rank_bundle(exterior_power(T_B3, 3)) == 10
    @test rank_bundle(exterior_power(T_B3, 4)) == 5
    @test rank_bundle(exterior_power(T_B3, 5)) == 1

    # D4/P1, dim=6: ranks = binomial(6, k)
    T_D4 = tangent_bundle(partial_flag_variety(TypeD{4}, (1,)))
    @test rank_bundle(exterior_power(T_D4, 1)) == 6
    @test rank_bundle(exterior_power(T_D4, 2)) == 15
    @test rank_bundle(exterior_power(T_D4, 3)) == 20
    @test rank_bundle(exterior_power(T_D4, 6)) == 1

    # G2/P1, dim=5: ranks = binomial(5, k)
    T_G2 = tangent_bundle(partial_flag_variety(TypeG2, (1,)))
    @test rank_bundle(exterior_power(T_G2, 1)) == 5
    @test rank_bundle(exterior_power(T_G2, 2)) == 10
    @test rank_bundle(exterior_power(T_G2, 5)) == 1

    # C3/P1, dim=5: ranks = binomial(5, k)
    T_C3 = tangent_bundle(partial_flag_variety(TypeC{3}, (1,)))
    @test rank_bundle(exterior_power(T_C3, 1)) == 5
    @test rank_bundle(exterior_power(T_C3, 2)) == 10
    @test rank_bundle(exterior_power(T_C3, 5)) == 1
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  H⁰ of exterior powers of tangent bundles (from reference output files)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "H⁰(∧ᵏT): A-types" begin
    # A3/P2 = Gr(2,4), reference: A3-P2.txt
    T = tangent_bundle(Gr(2, 4))
    @test dimensions(exterior_power(T, 0))[0] == 1
    @test dimensions(exterior_power(T, 1))[0] == 15
    @test dimensions(exterior_power(T, 2))[0] == 90
    @test dimensions(exterior_power(T, 3))[0] == 175
    @test dimensions(exterior_power(T, 4))[0] == 105

    # A4/P2 = Gr(2,5), reference: A4-P2.txt
    T2 = tangent_bundle(Gr(2, 5))
    @test dimensions(exterior_power(T2, 1))[0] == 24
    @test dimensions(exterior_power(T2, 2))[0] == 252
    @test dimensions(exterior_power(T2, 3))[0] == 1248
    @test dimensions(exterior_power(T2, 6))[0] == 1176

    # A7/P3 = Gr(3,8), reference: A7-P3.txt
    T3 = tangent_bundle(Gr(3, 8))
    @test dimensions(exterior_power(T3, 1))[0] == 63
    @test dimensions(exterior_power(T3, 2))[0] == 1890
  end

  @testset "H⁰(∧ᵏT): B-types" begin
    # B3/P1, reference: B3-P1.txt
    T = tangent_bundle(partial_flag_variety(TypeB{3}, (1,)))
    @test dimensions(exterior_power(T, 1))[0] == 21
    @test dimensions(exterior_power(T, 2))[0] == 189
    @test dimensions(exterior_power(T, 3))[0] == 616
    @test dimensions(exterior_power(T, 4))[0] == 819
    @test dimensions(exterior_power(T, 5))[0] == 378
  end

  @testset "H⁰(∧ᵏT): C-types" begin
    # C3/P1, reference: C3-P1.txt
    T = tangent_bundle(partial_flag_variety(TypeC{3}, (1,)))
    @test dimensions(exterior_power(T, 1))[0] == 35
    @test dimensions(exterior_power(T, 2))[0] == 280
    @test dimensions(exterior_power(T, 3))[0] == 840
    @test dimensions(exterior_power(T, 4))[0] == 1050
    @test dimensions(exterior_power(T, 5))[0] == 462
  end

  @testset "H⁰(∧ᵏT): D-types" begin
    # D4/P1, reference: D4-P1.txt
    T = tangent_bundle(partial_flag_variety(TypeD{4}, (1,)))
    @test dimensions(exterior_power(T, 1))[0] == 28
    @test dimensions(exterior_power(T, 2))[0] == 350
    @test dimensions(exterior_power(T, 3))[0] == 1680
    @test dimensions(exterior_power(T, 6))[0] == 1386
  end

  @testset "H⁰(∧ᵏT): G2" begin
    # G2/P1, reference: G2-P1.txt
    T = tangent_bundle(partial_flag_variety(TypeG2, (1,)))
    @test dimensions(exterior_power(T, 1))[0] == 21
    @test dimensions(exterior_power(T, 2))[0] == 189
    @test dimensions(exterior_power(T, 3))[0] == 616
    @test dimensions(exterior_power(T, 4))[0] == 819
    @test dimensions(exterior_power(T, 5))[0] == 378
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  String constructors
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "String constructors" begin
    X1 = partial_flag_variety("A4", 1)
    @test dimension(X1) == 4
    @test picard_rank(X1) == 1

    X2 = partial_flag_variety("A3", [1, 3])
    @test picard_rank(X2) == 2

    X3 = partial_flag_variety("B3", 1)
    @test dimension(X3) == 5

    X4 = partial_flag_variety("D4", [1, 3])
    @test picard_rank(X4) == 2
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Integer widening: BigInt, Int32, etc.
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Integer widening" begin
    # BigInt nodes (vector form accepts <:Integer)
    X_big = partial_flag_variety(TypeA{4}, [BigInt(2)])
    @test dimension(X_big) == 6

    # Int32 degree for line_bundle
    X = projective_space(3)
    L = line_bundle(X, Int32(2))
    @test rank_bundle(L) == 1

    # Vector of BigInt for line_bundle
    X2 = partial_flag_variety(TypeA{3}, (1, 3))
    L2 = line_bundle(X2, BigInt[1, 2])
    @test rank_bundle(L2) == 1
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Bundle type hierarchy
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Bundle type hierarchy" begin
    X = Gr(2, 4)
    T = tangent_bundle(X)
    @test T isa Bundle
    @test T isa CompletelyReducibleBundle
    @test T isa Bundle{marked_type(X)}
    @test structure_sheaf(X) isa Bundle
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  FilteredBundle
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "FilteredBundle" begin
    # Grassmannian: tangent bundle is irreducible → 1 filtration step
    X = Gr(2, 4)
    F = filtered_tangent_bundle(X)
    @test F isa FilteredBundle
    @test F isa Bundle
    @test n_filtration_steps(F) == 1
    @test rank_bundle(F) == dimension(X)
    @test rank_bundle(total_bundle(F)) == dimension(X)

    # Projective space: tangent is irreducible → 1 step
    X2 = projective_space(3)
    F2 = filtered_tangent_bundle(X2)
    @test n_filtration_steps(F2) == 1
    @test rank_bundle(F2) == 3

    # Full flag: more steps
    X3 = full_flag_variety(TypeA{3})
    F3 = filtered_tangent_bundle(X3)
    @test n_filtration_steps(F3) >= 2
    @test rank_bundle(F3) == dimension(X3)

    # Total bundle of filtered tangent = tangent bundle (same components)
    X4 = Gr(2, 5)
    F4 = filtered_tangent_bundle(X4)
    T4 = tangent_bundle(X4)
    @test rank_bundle(total_bundle(F4)) == rank_bundle(T4)

    # Tensor product with structure sheaf preserves filtration
    S = structure_sheaf(X)
    FS = tensor_product(F, S)
    @test n_filtration_steps(FS) == n_filtration_steps(F)
    @test rank_bundle(FS) == rank_bundle(F)

    # Exterior powers of filtered bundles
    d = dimension(X)
    @test rank_bundle(exterior_power(F, 0)) == 1
    @test rank_bundle(exterior_power(F, 1)) == d
    @test rank_bundle(exterior_power(F, 2)) == binomial(d, 2)
    @test rank_bundle(exterior_power(F, d)) == 1

    # Symmetric powers of filtered bundles
    @test rank_bundle(symmetric_power(F, 0)) == 1
    @test rank_bundle(symmetric_power(F, 1)) == d
    @test rank_bundle(symmetric_power(F, 2)) == binomial(d + 1, 2)

    # Multi-step filtration: exterior power increases steps
    F3w2 = exterior_power(F3, 2)
    d3 = dimension(X3)
    @test rank_bundle(F3w2) == binomial(d3, 2)
    @test n_filtration_steps(F3w2) >= n_filtration_steps(F3)

    # Dual reverses filtration
    Fd = dual(F3)
    @test n_filtration_steps(Fd) == n_filtration_steps(F3)
    @test rank_bundle(Fd) == rank_bundle(F3)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Universal bundles
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Universal bundles" begin
    # Gr(2, 5): U has rank 2, Q has rank 3
    X = Gr(2, 5)
    U = universal_subbundle(X)
    Q = universal_quotient_bundle(X)
    @test rank_bundle(U) == 2
    @test rank_bundle(Q) == 3

    # Gr(3, 7)
    X2 = Gr(3, 7)
    @test rank_bundle(universal_subbundle(X2)) == 3
    @test rank_bundle(universal_quotient_bundle(X2)) == 4

    # ℙ^4 = Gr(1, 5)
    X3 = projective_space(4)
    @test rank_bundle(universal_subbundle(X3)) == 1
    @test rank_bundle(universal_quotient_bundle(X3)) == 4

    # Tautological bundles on partial flag
    # On Fl(1,2;4) = A₃/P_{1,2}, the irreducible equivariant bundles
    # E_{ω₁} and E_{ω₂} have fiber dimensions 1 and 1 respectively.
    # (These are NOT the geometric tautological subbundles V₁, V₂;
    #  the latter are filtered, not completely reducible.)
    X4 = flag_variety(4, (1, 2))
    Us = tautological_bundles(X4)
    @test length(Us) == 2
    @test rank_bundle(Us[1]) == 1
    @test rank_bundle(Us[2]) == 1
  end

  @testset "Spinor bundles" begin
    # Odd quadric Q^5 = B_3/P_1: single spinor of rank 2^{3-1} = 4
    X5 = quadric(5)
    S5 = spinor_bundle(X5)
    @test rank_bundle(S5) == 4

    # Even quadric Q^4 = D_3/P_1: two half-spinors of rank 2^{3-2} = 2
    X4 = quadric(4)
    Sp = spinor_bundle(X4, :plus)
    Sm = spinor_bundle(X4, :minus)
    @test rank_bundle(Sp) == 2
    @test rank_bundle(Sm) == 2

    # Both half-spinors together
    S4 = spinor_bundle(X4)
    @test rank_bundle(S4) == 4

    # Q^3 = B_2/P_1: single spinor of rank 2
    X3 = quadric(3)
    S3 = spinor_bundle(X3)
    @test rank_bundle(S3) == 2
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Hodge numbers
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Hodge numbers" begin
    # ℙ^3: diagonal Hodge diamond
    X = projective_space(3)
    H = hodge_numbers(X)
    @test H[1, 1] == 1  # h^{0,0}
    @test H[2, 2] == 1  # h^{1,1}
    @test H[3, 3] == 1  # h^{2,2}
    @test H[4, 4] == 1  # h^{3,3}
    @test H[1, 2] == 0  # h^{0,1}
    @test H[2, 1] == 0  # h^{1,0}

    # Gr(2, 4): b₀=1, b₂=1, b₄=2, b₆=1, b₈=1
    #           h^{p,q} = 0 for p ≠ q, h^{p,p} = b_{2p}
    X2 = Gr(2, 4)
    H2 = hodge_numbers(X2)
    @test H2[1, 1] == 1  # h^{0,0} = b₀ = 1
    @test H2[2, 2] == 1  # h^{1,1} = b₂ = 1
    @test H2[3, 3] == 2  # h^{2,2} = b₄ = 2
    @test H2[1, 3] == 0
  end

  @testset "Twisted Hodge numbers" begin
    # ℙ³, twist 0: h^q(Ω^p) for p ≠ q is 0
    X = projective_space(3)
    H0 = twisted_hodge_numbers(X, 0)
    @test H0[1, 1] == 1  # h^0(𝒪)
    @test H0[2, 2] == 1  # h^1(Ω^1)
    @test H0[3, 3] == 1  # h^2(Ω^2)
    @test H0[4, 4] == 1  # h^3(Ω^3)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Hochschild cohomology
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Hochschild cohomology" begin
    # ℙ^2: HH^0 = 1, HH^2 = h^0(∧²T) + h^1(T) + h^2(𝒪)
    X = projective_space(2)
    P = hochschild_cohomology(X)
    @test P isa PolyvectorParallelogram
    @test P[0, 0] == 1    # h^0(𝒪) = 1
    @test P[1, 0] == 8    # h^0(T_ℙ²) = h^0(𝒪(1)^3) = dim Aut(ℙ²) = 8
    @test P[2, 0] == 10   # h^0(∧²T) = h^0(𝒪(3)) = 10
    @test P[0, 1] == 0    # h^1(𝒪) = 0
    @test P[0, 2] == 0    # h^2(𝒪) = 0

    # ℙ^1
    X1 = projective_space(1)
    P1 = hochschild_cohomology(X1)
    @test P1[0, 0] == 1   # h^0(𝒪) = 1
    @test P1[1, 0] == 3   # h^0(T_ℙ¹) = dim SL₂ = 3
    @test P1[0, 1] == 0   # h^1(𝒪) = 0
    @test P1[1, 1] == 0   # h^1(T_ℙ¹) = 0
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Koszul: SES solver
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Koszul SES solver" begin
    # Simple case: 0 → A → B → C → 0 on dim=1 variety
    # H*(A) = [1,0], H*(B) = [2,0]
    # H*(C) should be [1,0] (determined)
    a = Cohomology{BigInt}(BigInt[1, 0], 1)
    b = Cohomology{BigInt}(BigInt[2, 0], 1)
    (c, det) = solve_ses_cohomology(a, b)
    @test c[0] == 1
    @test c[1] == 0
    @test det == true

    # Case with connecting homomorphism ambiguity
    # H*(A) = [0,0,1], H*(B) = [0,0,0] on dim=2 variety
    a2 = Cohomology{BigInt}(BigInt[0, 0, 1], 2)
    b2 = Cohomology{BigInt}(BigInt[0, 0, 0], 2)
    (c2, det2) = solve_ses_cohomology(a2, b2)
    # χ(C) = χ(B) - χ(A) = 0 - 1 = -1
    @test sum((-1)^i * c2[i] for i in 0:2) == -1
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Zero loci
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "ZeroLocus: construction" begin
    X = projective_space(4)
    E = line_bundle(X, 5)
    Z = zero_locus(E)

    @test dimension(Z) == 3
    @test codimension(Z) == 1
    @test ambient_variety(Z) === X
  end

  @testset "ZeroLocus: Euler characteristic" begin
    # Quintic CY3: χ(O_Z) = 0
    X = projective_space(4)
    Z = zero_locus(line_bundle(X, 5))
    @test euler_characteristic(Z) == 0

    # Quartic K3: χ(O_Z) = 2
    X3 = projective_space(3)
    Z3 = zero_locus(line_bundle(X3, 4))
    @test euler_characteristic(Z3) == 2

    # Quadric surface in P^3: χ(O_Z) = 1
    Z_q = zero_locus(line_bundle(X3, 2))
    @test euler_characteristic(Z_q) == 1
  end

  @testset "ZeroLocus: CY detection" begin
    X = projective_space(4)

    # Quintic in P^4: CY candidate and CY
    @test is_calabi_yau_candidate(line_bundle(X, 5)) == true
    @test is_calabi_yau(zero_locus(line_bundle(X, 5))) == true

    # Cubic in P^4: NOT CY candidate (degree 3 ≠ index 5)
    @test is_calabi_yau_candidate(line_bundle(X, 3)) == false

    # Two cubics in P^5: CY candidate and CY
    X5 = projective_space(5)
    E = direct_sum(line_bundle(X5, 3), line_bundle(X5, 3))
    @test is_calabi_yau_candidate(E) == true
    @test is_calabi_yau(zero_locus(E)) == true
  end

  @testset "ZeroLocus: Hodge numbers (CY3)" begin
    # Quintic CY3: h^{1,1}=1, h^{2,1}=101
    X = projective_space(4)
    Z = zero_locus(line_bundle(X, 5))
    h = hodge_numbers(Z)
    @test h[1, 1] == 1   # h^{0,0}
    @test h[1, 4] == 1   # h^{0,3}
    @test h[2, 2] == 1   # h^{1,1}
    @test h[3, 2] == 101 # h^{2,1}
    @test h[2, 3] == 101 # h^{1,2} = h^{2,1}
    @test h[3, 3] == 1   # h^{2,2} = h^{1,1}

    # Two cubics in P^5: CY3, h^{1,1}=1, h^{2,1}=73
    X5 = projective_space(5)
    E = direct_sum(line_bundle(X5, 3), line_bundle(X5, 3))
    Z5 = zero_locus(E)
    h5 = hodge_numbers(Z5)
    @test h5[2, 2] == 1   # h^{1,1}
    @test h5[3, 2] == 73  # h^{2,1}
  end

  @testset "ZeroLocus: Hodge numbers (Fano 4-fold)" begin
    X = projective_space(5)

    # Quadric 4-fold: h^{1,1}=1, h^{2,2}=2
    Z_q = zero_locus(line_bundle(X, 2))
    h_q = hodge_numbers(Z_q)
    @test h_q[2, 2] == 1  # h^{1,1}
    @test h_q[3, 3] == 2  # h^{2,2}

    # Cubic 4-fold: h^{1,1}=1, h^{1,3}=1, h^{2,2}=21
    Z_c = zero_locus(line_bundle(X, 3))
    h_c = hodge_numbers(Z_c)
    @test h_c[2, 2] == 1  # h^{1,1}
    @test h_c[2, 4] == 1  # h^{1,3}
    @test h_c[3, 3] == 21 # h^{2,2}

    # Quartic 4-fold: h^{1,1}=1, h^{1,3}=21, h^{2,2}=142
    Z_4 = zero_locus(line_bundle(X, 4))
    h_4 = hodge_numbers(Z_4)
    @test h_4[2, 2] == 1   # h^{1,1}
    @test h_4[2, 4] == 21  # h^{1,3}
    @test h_4[3, 3] == 142 # h^{2,2}
  end

  @testset "ZeroLocus: Hodge numbers (Grassmannian)" begin
    # O(1)^4 on Gr(2,6): Küchle c6, h^{1,1}=1, h^{2,2}=8
    X = Gr(2, 6)
    E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:4])
    Z = zero_locus(E)
    @test dimension(Z) == 4
    h = hodge_numbers(Z)
    @test h[2, 2] == 1  # h^{1,1}
    @test h[3, 3] == 8  # h^{2,2}
  end

  @testset "ZeroLocus: conormal χ recursion" begin
    # Verify χ(Ω^p_Z) via the conormal recursion agrees with
    # the known identity χ(Ω^p) = (-1)^d χ(Ω^{d-p}) (Serre)
    X = projective_space(5)
    Z = zero_locus(line_bundle(X, 2))
    d = dimension(Z)  # 4

    chi = [PartialFlagVarieties._chi_omega_p_conormal(Z, p) for p in 0:d]
    # χ(Ω^0) = 1, χ(Ω^4) = 1 (Serre)
    @test chi[1] == chi[5]
    # χ(Ω^1) = χ(Ω^3) (Serre, d=4, (-1)^4=1)
    @test chi[2] == chi[4]
  end

  @testset "ZeroLocus: quartic K3 Hodge" begin
    X = projective_space(3)
    Z = zero_locus(line_bundle(X, 4))
    @test dimension(Z) == 2
    h = hodge_numbers(Z)
    @test h[1, 1] == 1   # h^{0,0}
    @test h[2, 2] == 20  # h^{1,1}
    @test h[1, 3] == 1   # h^{0,2}
    @test h[3, 1] == 1   # h^{2,0}
    @test h[3, 3] == 1   # h^{2,2} = h^{0,0}
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Exceptional collections
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "is_exceptional: basic bundles" begin
    X = projective_space(3)
    @test is_exceptional(structure_sheaf(X))
    @test is_exceptional(line_bundle(X, 1))
    @test is_exceptional(line_bundle(X, -1))
    @test is_exceptional(tangent_bundle(X))
    @test is_exceptional(cotangent_bundle(X))
  end

  @testset "is_exceptional_pair: line bundles on ℙ³" begin
    X = projective_space(3)
    O0 = structure_sheaf(X)
    O1 = line_bundle(X, 1)
    O2 = line_bundle(X, 2)
    # (Oₐ, O_b) is exceptional pair iff a < b
    @test is_exceptional_pair(O0, O1)
    @test is_exceptional_pair(O0, O2)
    @test is_exceptional_pair(O1, O2)
    @test !is_exceptional_pair(O1, O0)
  end

  @testset "is_strong_exceptional_pair: line bundles on ℙ³" begin
    X = projective_space(3)
    O0 = structure_sheaf(X)
    O1 = line_bundle(X, 1)
    @test is_strong_exceptional_pair(O0, O1)
  end

  @testset "Beilinson collection on ℙ⁴" begin
    X = projective_space(4)
    Es = beilinson_collection(X)
    @test length(Es) == 5
    @test is_full_exceptional_sequence(Es, X)
    @test is_strong_exceptional_sequence(Es)
    # dual Beilinson
    Ed = beilinson_collection_dual(X)
    @test length(Ed) == 5
    @test is_full_exceptional_sequence(Ed, X)
    @test is_strong_exceptional_sequence(Ed)
  end

  @testset "Kapranov collection on quadrics" begin
    # Odd quadrics Q^{2m-1}: B_m/P_1, collection has n+1 elements
    Q3 = quadric(3)   # B_2/P_1
    Es3 = kapranov_collection(Q3)
    @test length(Es3) == 4   # = chi(Q^3)
    @test is_full_exceptional_sequence(Es3, Q3)

    Q5 = quadric(5)   # B_3/P_1
    Es5 = kapranov_collection(Q5)
    @test length(Es5) == 6   # = chi(Q^5)
    @test is_full_exceptional_sequence(Es5, Q5)

    # Even quadrics Q^{2m}: D_m/P_1, collection has n+2 elements
    Q4 = quadric(4)   # D_3/P_1
    Es4 = kapranov_collection(Q4)
    @test length(Es4) == 6   # = chi(Q^4)
    @test is_full_exceptional_sequence(Es4, Q4)
  end

  @testset "schur_functor ranks on Gr(2,4)" begin
    X = Gr(2, 4)
    @test rank_bundle(schur_functor(X, [0, 0])) == 1     # O
    @test rank_bundle(schur_functor(X, [1, 0])) == 2     # U^∨
    @test rank_bundle(schur_functor(X, [1, 1])) == 1     # det(U^∨)
    @test rank_bundle(schur_functor(X, [2, 0])) == 3     # Sym²(U^∨)
    @test rank_bundle(schur_functor(X, [2, 1])) == 2     # U^∨ ⊗ det
    @test rank_bundle(schur_functor(X, [2, 2])) == 1     # det²
    # NTuple form
    @test rank_bundle(schur_functor(X, (2, 0))) == 3
  end

  @testset "schur_functor ranks on Gr(2,5)" begin
    X = Gr(2, 5)
    # GL(2) Schur functors: Sym^d has rank d+1
    @test rank_bundle(schur_functor(X, [0, 0])) == 1
    @test rank_bundle(schur_functor(X, [1, 0])) == 2
    @test rank_bundle(schur_functor(X, [1, 1])) == 1
    @test rank_bundle(schur_functor(X, [2, 0])) == 3
    @test rank_bundle(schur_functor(X, [3, 0])) == 4
    @test rank_bundle(schur_functor(X, [3, 3])) == 1
  end

  @testset "Kapranov collection on Gr(2,4)" begin
    X = Gr(2, 4)
    Es = kapranov_bundles_grassmannian(X)
    @test length(Es) == euler_characteristic(X)
    @test is_full_exceptional_sequence(Es, X)
    @test is_strong_exceptional_sequence(Es)
  end

  @testset "Kapranov collection on Gr(2,5)" begin
    X = Gr(2, 5)
    Es = kapranov_bundles_grassmannian(X)
    @test length(Es) == euler_characteristic(X)   # = 10
    @test is_full_exceptional_sequence(Es, X)
    @test is_strong_exceptional_sequence(Es)
  end

  @testset "Kapranov collection on Gr(3,6)" begin
    X = Gr(3, 6)
    Es = kapranov_bundles_grassmannian(X)
    @test length(Es) == euler_characteristic(X)   # = 20
    @test is_full_exceptional_sequence(Es, X)
    @test is_strong_exceptional_sequence(Es)
  end

end  # @testset "PartialFlagVarieties.jl"
