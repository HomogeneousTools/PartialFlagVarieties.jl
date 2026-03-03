using Test
using PartialFlagVarieties
using Lie
using StaticArrays

@testset "PartialFlagVarieties.jl" begin

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
  #  Dimensions
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Dimensions" begin
    # Projective spaces: dim(ℙⁿ) = n
    @test dimension(MarkedDynkinType{TypeA{1},(1,)}) == 1
    @test dimension(MarkedDynkinType{TypeA{2},(1,)}) == 2
    @test dimension(MarkedDynkinType{TypeA{4},(1,)}) == 4

    # Grassmannians: dim(Gr(k,n)) = k(n-k)
    @test dimension(MarkedDynkinType{TypeA{3},(2,)}) == 4    # Gr(2,4)
    @test dimension(MarkedDynkinType{TypeA{4},(2,)}) == 6    # Gr(2,5)
    @test dimension(MarkedDynkinType{TypeA{5},(3,)}) == 9    # Gr(3,6)

    # Full flag: dim(G/B) = |Φ⁺|
    @test dimension(MarkedDynkinType{TypeA{2},(1,2)}) == 3   # A2: 3 pos roots
    @test dimension(MarkedDynkinType{TypeA{3},(1,2,3)}) == 6 # A3: 6 pos roots

    # Exceptional: Cayley plane OP² = E6/P1 has dim 16
    @test dimension(MarkedDynkinType{TypeE{6},(1,)}) == 16

    # Quadrics: Q_n has dim n
    @test dimension(MarkedDynkinType{TypeB{2},(1,)}) == 3    # Q3
    @test dimension(MarkedDynkinType{TypeD{3},(1,)}) == 4    # Q4

    # OGr(5,10) = D5/P5: spinor variety, dim 10
    @test dimension(MarkedDynkinType{TypeD{5},(5,)}) == 10
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Euler characteristics
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Euler characteristics" begin
    # χ(ℙⁿ) = n + 1
    @test euler_characteristic(MarkedDynkinType{TypeA{4},(1,)}) == 5

    # χ(Gr(k,n)) = C(n,k)
    @test euler_characteristic(MarkedDynkinType{TypeA{3},(2,)}) == 6    # C(4,2)
    @test euler_characteristic(MarkedDynkinType{TypeA{4},(2,)}) == 10   # C(5,2)
    @test euler_characteristic(MarkedDynkinType{TypeA{5},(2,)}) == 15   # C(6,2)
    @test euler_characteristic(MarkedDynkinType{TypeA{5},(3,)}) == 20   # C(6,3)

    # χ(OP²) = 27 (Cayley plane)
    @test euler_characteristic(MarkedDynkinType{TypeE{6},(1,)}) == 27

    # χ(E7/P7) = 56 (Freudenthal variety)
    @test euler_characteristic(MarkedDynkinType{TypeE{7},(7,)}) == 56

    # χ(OGr(5,10)) = 2⁴ = 16
    @test euler_characteristic(MarkedDynkinType{TypeD{5},(5,)}) == 16

    # χ(G/B) = |W|
    @test euler_characteristic(MarkedDynkinType{TypeA{2},(1,2)}) == 6
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Betti numbers
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Betti numbers" begin
    # ℙⁿ: all Betti numbers = 1
    betti_P2 = betti_numbers(MarkedDynkinType{TypeA{2},(1,)})
    @test betti_P2 == [1, 1, 1]

    betti_P4 = betti_numbers(MarkedDynkinType{TypeA{4},(1,)})
    @test betti_P4 == [1, 1, 1, 1, 1]

    # Gr(2,4): Betti = 1, 1, 2, 1, 1
    betti_Gr24 = betti_numbers(MarkedDynkinType{TypeA{3},(2,)})
    @test betti_Gr24 == [1, 1, 2, 1, 1]

    # Gr(2,5): Betti = 1, 1, 2, 2, 2, 1, 1
    betti_Gr25 = betti_numbers(MarkedDynkinType{TypeA{4},(2,)})
    @test betti_Gr25 == [1, 1, 2, 2, 2, 1, 1]

    # Sum of Betti = χ
    betti_OP2 = betti_numbers(MarkedDynkinType{TypeE{6},(1,)})
    @test sum(betti_OP2) == 27
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Special matrix
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Special matrix" begin
    MDT = MarkedDynkinType{TypeA{3},(2,)}
    M = special_matrix(MDT)
    Minv = special_matrix_inv(MDT)

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
    @test is_generalized_grassmannian(MarkedDynkinType{TypeA{3},(2,)})
    @test !is_generalized_grassmannian(MarkedDynkinType{TypeA{3},(1,3)})

    @test is_full_flag(MarkedDynkinType{TypeA{2},(1,2)})
    @test !is_full_flag(MarkedDynkinType{TypeA{3},(2,)})

    # Cominuscule
    @test is_cominuscule(MarkedDynkinType{TypeA{4},(2,)})      # Gr(2,5)
    @test is_cominuscule(MarkedDynkinType{TypeD{5},(1,)})      # Quadric
    @test is_cominuscule(MarkedDynkinType{TypeD{5},(5,)})      # Spinor
    @test is_cominuscule(MarkedDynkinType{TypeE{6},(1,)})      # OP²
    @test !is_cominuscule(MarkedDynkinType{TypeB{3},(2,)})

    # Minuscule
    @test is_minuscule(MarkedDynkinType{TypeA{4},(2,)})
    @test is_minuscule(MarkedDynkinType{TypeB{3},(3,)})
    @test is_minuscule(MarkedDynkinType{TypeD{5},(5,)})
    @test !is_minuscule(MarkedDynkinType{TypeB{3},(1,)})

    # Adjoint
    @test is_adjoint(MarkedDynkinType{TypeA{3},(1,3)})         # ℙ(T*ℙ³)
    @test is_adjoint(MarkedDynkinType{TypeB{3},(2,)})
    @test is_adjoint(MarkedDynkinType{TypeG2,(2,)})
    @test is_adjoint(MarkedDynkinType{TypeE{6},(2,)})

    # Coadjoint
    @test is_coadjoint(MarkedDynkinType{TypeB{3},(1,)})
    @test is_coadjoint(MarkedDynkinType{TypeG2,(1,)})
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
    @test length(nonpar) == dimension(MDT)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  PartialFlagVariety wrapper
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "PartialFlagVariety wrapper" begin
    V = partial_flag_variety(TypeA{3}, (2,))
    @test dimension(V) == 4
    @test euler_characteristic(V) == 6
    @test marked_nodes(V) == (2,)

    V_full = full_flag_variety(TypeA{2})
    @test is_full_flag(V_full)
    @test dimension(V_full) == 3
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
    @test d.central == -rep.central
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  CompletelyReducibleBundle
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Structure sheaf" begin
    MDT = MarkedDynkinType{TypeA{3},(2,)}
    O = structure_sheaf(MDT)
    @test rank_bundle(O) == 1
    @test n_components(O) == 1
  end

  @testset "Line bundles" begin
    MDT = MarkedDynkinType{TypeA{4},(1,)}
    L = line_bundle(MDT, 1)
    @test rank_bundle(L) == 1
  end

  @testset "Tangent and cotangent bundles" begin
    # Gr(2,4): T = S* ⊗ Q, both have rank 2, so T has rank 4
    MDT = MarkedDynkinType{TypeA{3},(2,)}
    T = tangent_bundle(MDT)
    Ω = cotangent_bundle(MDT)
    @test rank_bundle(T) == 4
    @test rank_bundle(Ω) == 4

    # ℙ⁴: T has rank 4
    MDT2 = MarkedDynkinType{TypeA{4},(1,)}
    T2 = tangent_bundle(MDT2)
    @test rank_bundle(T2) == 4
  end

  @testset "Bundle operations" begin
    MDT = MarkedDynkinType{TypeA{4},(1,)}
    O = structure_sheaf(MDT)
    T = tangent_bundle(MDT)

    # O ⊗ T = T
    @test rank_bundle(tensor_product(O, T)) == rank_bundle(T)

    # ⊕
    @test rank_bundle(direct_sum(T, O)) == rank_bundle(T) + 1

    # ⋀⁰T = O, ⋀¹T = T
    @test rank_bundle(exterior_power(T, 0)) == 1
    @test rank_bundle(exterior_power(T, 1)) == rank_bundle(T)
  end

  @testset "Twist" begin
    MDT = MarkedDynkinType{TypeA{4},(1,)}
    O = structure_sheaf(MDT)
    Ot = twist(O, 1, 3)
    @test rank_bundle(Ot) == 1
  end

  @testset "Determinant bundle" begin
    MDT = MarkedDynkinType{TypeA{3},(2,)}
    T = tangent_bundle(MDT)
    det_T = det_bundle(T)
    @test rank_bundle(det_T) == 1
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Cohomology
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Cohomology: H*(ℙⁿ, 𝒪)" begin
    MDT = MarkedDynkinType{TypeA{4},(1,)}
    O = structure_sheaf(MDT)
    H = cohomology(MDT, O)
    d = dimensions(H)

    # H⁰(ℙ⁴, 𝒪) = 1, Hⁱ = 0 for i > 0
    @test d[0] == 1
    for i in 1:4
      @test d[i] == 0
    end
  end

  @testset "Cohomology: H*(ℙⁿ, 𝒪(k))" begin
    MDT = MarkedDynkinType{TypeA{4},(1,)}

    # H⁰(ℙ⁴, 𝒪(1)) = 5 (standard rep of SL(5))
    L1 = line_bundle(MDT, 1)
    H1 = dimensions(cohomology(MDT, L1))
    @test H1[0] == 5

    # H⁰(ℙ⁴, 𝒪(2)) = C(6,2) = 15
    L2 = twist(structure_sheaf(MDT), 1, 2)
    H2 = dimensions(cohomology(MDT, L2))
    @test H2[0] == 15

    # H⁰(ℙ⁴, 𝒪(3)) = C(7,3) = 35
    L3 = twist(structure_sheaf(MDT), 1, 3)
    H3 = dimensions(cohomology(MDT, L3))
    @test H3[0] == 35
  end

  @testset "Cohomology: Euler characteristic" begin
    MDT = MarkedDynkinType{TypeA{4},(1,)}

    # χ(ℙ⁴, 𝒪) = 1
    @test euler_char_bundle(MDT, structure_sheaf(MDT)) == 1

    # χ(ℙ⁴, 𝒪(1)) = 5
    @test euler_char_bundle(MDT, line_bundle(MDT, 1)) == 5
  end

  @testset "Cohomology: 0-based indexing" begin
    MDT = MarkedDynkinType{TypeA{2},(1,)}
    O = structure_sheaf(MDT)
    H = dimensions(cohomology(MDT, O))
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
      @test dimension(MDT) == length(positive_nonparabolic_roots(MDT))
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
      MDT = MarkedDynkinType{DT,marks}
      @test euler_characteristic(MDT) == sum(betti_numbers(MDT))
    end
  end

  @testset "Consistency: rank(T) = dim(X)" begin
    for (DT, marks) in [
      (TypeA{3}, (2,)),
      (TypeA{4}, (1,)),
      (TypeB{3}, (1,)),
    ]
      MDT = MarkedDynkinType{DT,marks}
      @test Int(rank_bundle(tangent_bundle(MDT))) == dimension(MDT)
    end
  end

end  # @testset "PartialFlagVarieties.jl"
