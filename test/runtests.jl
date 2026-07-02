using Test
using PartialFlagVarieties
using Semisimple
using StaticArrays

mdt(::Type{DT}, marked) where {DT<:DynkinType} = MarkedDynkinType(DT, marked)

@testset "PartialFlagVarieties.jl" begin

  # ═══════════════════════════════════════════════════════════════════════════
  #  Semisimple.jl extensions: cartan_type, parse_dynkin_type
  #
  #  REMINDER: cartan_type, cartan_type_with_ordering, and parse_dynkin_type are
  #  defined in src/Semisimple.jl as staging code "for inclusion in the
  #  Semisimple.jl package" (see that file's header). They are currently owned
  #  and exported by PartialFlagVarieties, so their tests live here. When these
  #  functions are upstreamed into Semisimple.jl, the cartan_type and
  #  parse_dynkin_type testsets below should move into Semisimple.jl's own suite.
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
    @test DT === ProductDynkinType{Tuple{TypeA{2},TypeB{3}}}

    # Whitespace tolerance
    @test parse_dynkin_type(" A3 ") === TypeA{3}
  end

  @testset "parse_dynkin_type: malformed input" begin
    # Empty or whitespace-only input.
    @test_throws ArgumentError parse_dynkin_type("")
    @test_throws ArgumentError parse_dynkin_type("   ")
    # No rank digits, or digits before the letter.
    @test_throws ArgumentError parse_dynkin_type("A")
    @test_throws ArgumentError parse_dynkin_type("3A")
    # Letter outside the A–G family alphabet.
    @test_throws ArgumentError parse_dynkin_type("Z9")
    # Separator with no parseable component on either side.
    @test_throws ArgumentError parse_dynkin_type("xx")
    # Well-formed letter+digit, but not a valid Cartan rank for that family.
    @test_throws ArgumentError parse_dynkin_type("A0")  # A needs rank >= 1
    @test_throws ArgumentError parse_dynkin_type("B1")  # B needs rank >= 2
    @test_throws ArgumentError parse_dynkin_type("D3")  # D needs rank >= 4
    @test_throws ArgumentError parse_dynkin_type("E5")  # E only ranks 6, 7, 8
    @test_throws ArgumentError parse_dynkin_type("F3")  # F only rank 4
    @test_throws ArgumentError parse_dynkin_type("G3")  # G only rank 2
    # The same guard fires when a bad string reaches a variety constructor.
    @test_throws ArgumentError partial_flag_variety("Q7", 1)
    @test_throws ArgumentError partial_flag_variety("", 1)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  MarkedDynkinType
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "MarkedDynkinType basics" begin
    MDT = mdt(TypeA{4}, (2,))

    @test marked_nodes(MDT) == (2,)
    @test unmarked_nodes(MDT) == (1, 3, 4)
    @test central_rank(MDT) == 1
    @test levi_rank(MDT) == 3
    @test rank(MDT) == 4

    # Full flag
    MDT_full = mdt(TypeA{3}, (1, 2, 3))
    @test levi_type(MDT_full) === nothing
    @test central_rank(MDT_full) == 3
    @test levi_rank(MDT_full) == 0

    # Two marked nodes
    MDT2 = mdt(TypeA{4}, (1, 3))
    @test marked_nodes(MDT2) == (1, 3)
    @test unmarked_nodes(MDT2) == (2, 4)
    @test central_rank(MDT2) == 2
  end

  @testset "MarkedDynkinType constructors" begin
    mdt1 = MarkedDynkinType(TypeA{3}, (2,))
    @test dynkin_type(mdt1) === TypeA{3}
    @test marked_nodes(mdt1) == (2,)

    mdt2 = MarkedDynkinType(TypeB{4}, [1, 3])
    @test dynkin_type(mdt2) === TypeB{4}
    @test marked_nodes(mdt2) == (1, 3)

    mdt3 = MarkedDynkinType(TypeD{5}, 5)
    @test dynkin_type(mdt3) === TypeD{5}
    @test marked_nodes(mdt3) == (5,)
  end

  @testset "is_borel" begin
    @test is_borel(marked_dynkin_type(full_flag_variety(TypeA{3}))) == true
    @test is_borel(marked_dynkin_type(full_flag_variety(TypeB{3}))) == true
    @test is_borel(marked_dynkin_type(full_flag_variety(TypeG2))) == true
    @test is_borel(marked_dynkin_type(Gr(2, 4))) == false
    @test is_borel(marked_dynkin_type(projective_space(3))) == false
    @test is_borel(marked_dynkin_type(adjoint_variety(TypeG2))) == false
  end

  @testset "central_scaling_factor" begin
    @test central_scaling_factor(marked_dynkin_type(projective_space(3))) == 4
    @test central_scaling_factor(marked_dynkin_type(projective_space(4))) == 5
    @test central_scaling_factor(marked_dynkin_type(Gr(2, 4))) == 2

    @test central_scaling_factor(marked_dynkin_type(flag_variety(4, [1, 3]))) > 0
    @test central_scaling_factor(marked_dynkin_type(quadric(5))) > 0
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Levi type identification
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Levi type computation" begin
    # Gr(2,5) = A4/P2 → Levi = A1 × A2
    @test levi_type(mdt(TypeA{4}, (2,))) ==
      ProductDynkinType{Tuple{TypeA{1},TypeA{2}}}

    # Gr(1,5) = A4/P1 → Levi = A3
    @test levi_type(mdt(TypeA{4}, (1,))) == TypeA{3}

    # D5/P5 → Levi = A4
    @test levi_type(mdt(TypeD{5}, (5,))) == TypeA{4}

    # B3/P1 → Levi = B2
    @test levi_type(mdt(TypeB{3}, (1,))) == TypeB{2}

    # E6/P1 → Levi = D5
    @test levi_type(mdt(TypeE{6}, (1,))) == TypeD{5}

    # E6/P2 → Levi = A4 × ... (compute and check rank)
    @test levi_rank(mdt(TypeE{6}, (2,))) == 5
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  PartialFlagVariety constructors
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "PartialFlagVariety constructors" begin
    V = partial_flag_variety(TypeA{3}, (2,))
    @test V isa PartialFlagVariety
    @test marked_nodes(V) == (2,)
    @test marked_dynkin_type(V) == MarkedDynkinType(TypeA{3}, (2,))
    @test unmarked_nodes(V) == (1, 3)
    @test dynkin_type(V) == TypeA{3}
    @test levi_type(V) == ProductDynkinType{Tuple{TypeA{1},TypeA{1}}}
    @test levi_rank(V) == 2
    @test central_rank(V) == 1

    V_full = full_flag_variety(TypeA{2})
    @test is_full_flag_variety(V_full)
    @test levi_type(V_full) === nothing

    # String constructor
    V_str = PartialFlagVariety("A3", [2])
    @test dimension(V_str) == 4
    @test marked_nodes(V_str) == (2,)

    # Product type string constructor
    V_prod = PartialFlagVariety("A2xB3", [1, 4])
    @test V_prod isa PartialFlagVariety
  end

  @testset "flag_variety dimension validation" begin
    # Valid, sorted, distinct, in-range dimensions work.
    @test marked_nodes(flag_variety(4, [1, 2])) == (1, 2)
    # n - 1 is the largest admissible dimension.
    @test marked_nodes(flag_variety(5, [1, 2, 4])) == (1, 2, 4)
    # Unsorted, repeated, or out-of-range dimensions are rejected, not silently fixed.
    @test_throws ArgumentError flag_variety(4, [2, 1])
    @test_throws ArgumentError flag_variety(4, [1, 1])
    @test_throws ArgumentError flag_variety(4, [0, 2])
    @test_throws ArgumentError flag_variety(4, [1, 4])
  end

  @testset "Constructor argument validation" begin
    # Grassmannian Gr(k, n): need 1 <= k <= n - 1.
    @test_throws ArgumentError Gr(0, 5)   # k too small
    @test_throws ArgumentError Gr(5, 5)   # k = n
    @test_throws ArgumentError Gr(6, 5)   # k > n
    @test_throws ArgumentError Gr(2, 2)   # k > n - 1

    # Projective space needs n >= 1.
    @test_throws ArgumentError projective_space(0)
    @test_throws ArgumentError projective_space(-1)

    # Orthogonal Grassmannian OGr(k, n): k >= 1 and k <= floor(n / 2).
    @test_throws ArgumentError OGr(0, 7)    # k too small
    @test_throws ArgumentError OGr(4, 7)    # B_3: need k <= 3
    @test_throws ArgumentError OGr(6, 10)   # D_5: need k <= 5

    # Symplectic Grassmannian SGr(k, n): n even and 1 <= k <= n / 2.
    @test_throws ArgumentError SGr(2, 7)    # n odd
    @test_throws ArgumentError SGr(4, 6)    # C_3: need k <= 3

    # Quadric needs n >= 1.
    @test_throws ArgumentError quadric(0)

    # Marked nodes must be in range and distinct (validated, not silently fixed).
    @test_throws ArgumentError partial_flag_variety(TypeA{3}, [4])     # node > rank
    @test_throws ArgumentError partial_flag_variety(TypeA{3}, [0])     # node < 1
    @test_throws ArgumentError partial_flag_variety(TypeA{3}, [2, 2])  # repeated node
    @test_throws ArgumentError partial_flag_variety("B3", [4])         # node > rank
  end

  @testset "PartialFlagVariety products" begin
    X = projective_space(1)
    Y = projective_space(2)
    XY = product(X, Y)

    @test XY isa PartialFlagVariety
    @test XY == X * Y
    @test XY == partial_flag_variety(ProductDynkinType{Tuple{TypeA{1},TypeA{2}}}, (1, 2))
    @test dimension(XY) == dimension(X) + dimension(Y)
    @test picard_rank(XY) == picard_rank(X) + picard_rank(Y)
    @test marked_nodes(XY) == (1, 2)

    XYZ = product(X, X, X)
    @test XYZ == (X * X) * X
    @test dimension(XYZ) == 3
    @test picard_rank(XYZ) == 3
    @test marked_nodes(XYZ) == (1, 2, 3)
    @test zerolocus62_label(XYZ) == "111"
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

    # Low-dimensional quadrics are accidentally isomorphic to (products of) ℙⁿ:
    #   Q¹ = ℙ¹  and  Q² = ℙ¹ × ℙ¹
    @test quadric(1) == projective_space(1)
    @test quadric(2) == projective_space(1) * projective_space(1)

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
    MDT = mdt(TypeA{3}, (2,))
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

    @test is_full_flag_variety(full_flag_variety(TypeA{2}))
    @test !is_full_flag_variety(Gr(2, 4))

    # Projective space — every model, including accidental isomorphisms
    @test is_projective_space(projective_space(4))              # A₄/P₁
    @test is_projective_space(Gr(4, 5))                         # dual ℙ⁴ = A₄/P₄
    @test is_projective_space(SGr(1, 6))                        # C₃/P₁ = ℙ⁵
    @test is_projective_space(OGr(2, 5))                        # B₂/P₂ = ℙ³
    @test is_projective_space(OGr(2, 6))                        # D₃/P₂ = ℙ³
    @test is_projective_space(OGr(3, 6))                        # D₃/P₃ = ℙ³
    @test is_projective_space(quadric(1))                       # Q¹ = ℙ¹
    @test !is_projective_space(Gr(2, 5))                        # proper Grassmannian
    @test !is_projective_space(quadric(2))                      # ℙ¹ × ℙ¹, not a ℙⁿ
    @test !is_projective_space(quadric(4))                      # D₃/P₁ quadric ≠ ℙ³
    @test !is_projective_space(LGr(2, 4))                          # C₂/P₂ = Q³
    @test !is_projective_space(OGr(3, 7))                       # B₃/P₃ spinor variety

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
    @test length(tangent_weights(mdt(TypeA{4}, (1,)))) == 1

    # Gr(2,4) = A3/P2: T = S* ⊗ Q (irreducible under Levi) → 1 tangent weight
    @test length(tangent_weights(mdt(TypeA{3}, (2,)))) == 1

    # Full flag: each positive root is maximal → #tangent weights = dim
    @test length(tangent_weights(mdt(TypeA{2}, (1, 2)))) == 3
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Roots decomposition
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "Root decomposition" begin
    MDT = mdt(TypeA{3}, (2,))

    nonpar = positive_nonparabolic_roots(MDT)
    par = positive_parabolic_roots(MDT)

    # Total = n_positive_roots(A3) = 6
    @test length(nonpar) + length(par) == n_positive_roots(TypeA{3})

    # dim(G/P) = number of nonparabolic positive roots
    X = partial_flag_variety(TypeA{3}, (2,))
    @test length(nonpar) == dimension(X)
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  IrrepLevi
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "IrrepLevi construction" begin
    MDT = mdt(TypeA{3}, (2,))

    # Fundamental weight ω₁
    ω₁ = fundamental_weight(TypeA{3}, 1)
    rep = IrrepLevi(MDT, ω₁)
    @test length(central_part(rep)) == 1
    @test semisimple_part(rep) isa WeightLatticeElem

    # Round-trip: to_ambient → IrrepLevi → to_ambient
    λ_back = to_ambient_weight(MDT, rep)
    @test λ_back == ω₁
  end

  @testset "fiber_dimension" begin
    X = Gr(2, 5)
    S_rep = only(components(universal_subbundle(X)))
    @test fiber_dimension(S_rep) == 2

    Q_rep = only(components(universal_quotient_bundle(X)))
    @test fiber_dimension(Q_rep) == 3

    O_rep = only(components(structure_sheaf(X)))
    @test fiber_dimension(O_rep) == 1

    L_rep = only(components(line_bundle(X, 3)))
    @test fiber_dimension(L_rep) == 1

    X_flag = full_flag_variety(TypeA{2})
    L_full = only(components(line_bundle(X_flag, [1, 2])))
    @test fiber_dimension(L_full) == 1
  end

  @testset "IrrepLevi round-trip" begin
    MDT = mdt(TypeA{4}, (2,))

    for i in 1:4
      ω = fundamental_weight(TypeA{4}, i)
      rep = IrrepLevi(MDT, ω)
      @test to_ambient_weight(MDT, rep) == ω
    end
  end

  @testset "IrrepLevi dual" begin
    MDT = mdt(TypeA{3}, (2,))
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

  @testset "CompletelyReducibleBundle constructors" begin
    X = Gr(3, 6)
    λ = WeightLatticeElem(dynkin_type(X), [0, 1, -1, 0, 0])
    E_weight = CompletelyReducibleBundle(X, λ)
    E_weights = CompletelyReducibleBundle(X, [λ, WeightLatticeElem(dynkin_type(X))])
    E_components = CompletelyReducibleBundle(X, copy(components(E_weight)))
    E_coeffs = CompletelyReducibleBundle(X, [0, 1, -1, 0, 0])
    E_sum = CompletelyReducibleBundle(X, [[0, 1, -1, 0, 0], [0, 0, 0, 0, 0]])
    O = structure_sheaf(X)
    Y = Gr(2, 5)
    bad_component = IrrepLevi(marked_dynkin_type(Y), [1, 0, 0, 0])

    @test variety(E_weight) === X
    @test components(E_weight) == components(universal_subbundle(X))
    @test rank_bundle(E_weight) == 3
    @test components(E_weights) == components(E_sum)
    @test components(E_components) == components(E_weight)
    @test components(E_coeffs) == components(E_weight)
    @test components(E_sum) == vcat(components(E_weight), components(O))
    @test rank_bundle(E_sum) == 4
    @test_throws ArgumentError CompletelyReducibleBundle(X, [bad_component])
  end

  @testset "Mismatched ambient varieties" begin
    # Binary bundle operations require both operands on the same variety.
    E = structure_sheaf(projective_space(2))
    F = structure_sheaf(projective_space(3))
    @test_throws ArgumentError tensor_product(E, F)
    @test_throws ArgumentError direct_sum(E, F)

    # A filtered bundle's pieces must live on its own variety.
    @test_throws ArgumentError FilteredBundle(
      Gr(2, 4), [structure_sheaf(projective_space(3))]
    )
  end

  @testset "Structure sheaf and zero bundle" begin
    X = Gr(2, 4)
    O = structure_sheaf(X)
    O2 = PartialFlagVarieties.O(X)
    Z = zero_bundle(X)
    @test rank_bundle(O) == 1
    @test variety(O2) === X
    @test components(O) == components(O2)
    @test n_components(O) == 1
    @test variety(O) === X
    @test p_dominant_weight(only(components(O))) == WeightLatticeElem(dynkin_type(X))
    @test rank_bundle(Z) == 0
    @test isempty(components(Z))
    @test n_components(Z) == 0
    @test variety(Z) === X
    @test components(direct_sum(O, Z)) == components(O)
    @test components(direct_sum(Z, O)) == components(O)
  end

  @testset "Line bundles" begin
    X = projective_space(4)
    L = line_bundle(X, 1)
    @test rank_bundle(L) == 1
    @test variety(L) === X
    @test p_dominant_weight(only(components(L))) ==
      WeightLatticeElem(dynkin_type(X), [1, 0, 0, 0])
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
    @test p_dominant_weight(only(components(L2))) ==
      WeightLatticeElem(dynkin_type(X2), [2, 0, 1])

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

  @testset "Shorthand constructors" begin
    X = Gr(2, 5)

    # S / Q: subbundle and quotient (type A → CompletelyReducibleBundle, has ==)
    @test S(X) == universal_subbundle(X)
    @test Q(X) == universal_quotient_bundle(X)

    # T: tangent bundle
    @test T(X) == tangent_bundle(X)

    # O: the one-arg structure sheaf is unaffected by the new line-bundle methods
    @test O(X) == structure_sheaf(X)
    @test O(X, 2) == line_bundle(X, 2)

    # E: general equivariant bundle from a weight, coefficients, or a list of them
    λ = WeightLatticeElem(dynkin_type(X), [0, 1, 0, 0])
    @test E(X, λ) == CompletelyReducibleBundle(X, λ)
    @test E(X, [0, 1, 0, 0]) == CompletelyReducibleBundle(X, [0, 1, 0, 0])
    @test E(X, [[0, 1, 0, 0], [0, 0, 0, 0]]) ==
      CompletelyReducibleBundle(X, [[0, 1, 0, 0], [0, 0, 0, 0]])

    # multi-step flag Fl(1, 3; 4): indexed subbundle and per-node line-bundle degrees
    Xf = partial_flag_variety(TypeA{3}, (1, 3))
    @test S(Xf, 1) == universal_subbundle(Xf, 1)
    @test rank_bundle(S(Xf, 2)) == rank_bundle(universal_subbundle(Xf, 2))
    @test O(Xf, [2, 3]) == line_bundle(Xf, [2, 3])
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
    det_T = det(T)
    @test rank_bundle(det_T) == 1

    L = line_bundle(X, 3)
    @test det(L) == L                              # det of a line bundle is itself

    M = line_bundle(X, 2)
    @test det(direct_sum(L, M)) == tensor_product(L, M)  # det(L⊕M) = L⊗M
  end

  @testset "picard_degrees" begin
    @test picard_degrees(line_bundle(projective_space(4), 3)) == [3]
    @test picard_degrees(structure_sheaf(projective_space(4))) == [0]
    @test picard_degrees(anticanonical_bundle(Gr(2, 5))) == [5]
    @test picard_degrees(canonical_bundle(Gr(2, 5))) == [-5]

    X_flag = flag_variety(4, [1, 2])
    @test picard_degrees(line_bundle(X_flag, [1, -2])) == [1, -2]

    P3 = projective_space(3)
    @test_throws ArgumentError picard_degrees(
      direct_sum(line_bundle(P3, 1), line_bundle(P3, 2))
    )
    @test_throws ArgumentError picard_degrees(tangent_bundle(P3))
  end

  @testset "det(tangent_bundle) == anticanonical_bundle" begin
    det_tangent_is_anticanonical(X) = det(tangent_bundle(X)) == anticanonical_bundle(X)

    @test det_tangent_is_anticanonical(projective_space(1))
    @test det_tangent_is_anticanonical(projective_space(2))
    @test det_tangent_is_anticanonical(projective_space(4))
    @test det_tangent_is_anticanonical(Gr(2, 4))
    @test det_tangent_is_anticanonical(Gr(2, 5))
    @test det_tangent_is_anticanonical(Gr(2, 6))
    @test det_tangent_is_anticanonical(Gr(3, 6))
    @test det_tangent_is_anticanonical(Gr(3, 7))
    @test det_tangent_is_anticanonical(Gr(4, 8))
    @test det_tangent_is_anticanonical(flag_variety(4, [1, 2]))
    @test det_tangent_is_anticanonical(flag_variety(5, [1, 2]))
    @test det_tangent_is_anticanonical(flag_variety(5, [2, 3]))
    @test det_tangent_is_anticanonical(flag_variety(5, [1, 4]))
    @test det_tangent_is_anticanonical(flag_variety(6, [2, 4]))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeA{2}))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeA{3}))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeA{4}))
    @test det_tangent_is_anticanonical(quadric(3))
    @test det_tangent_is_anticanonical(quadric(5))
    @test det_tangent_is_anticanonical(quadric(7))
    @test det_tangent_is_anticanonical(OGr(2, 7))
    @test det_tangent_is_anticanonical(OGr(3, 7))
    @test det_tangent_is_anticanonical(OGr(2, 9))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeB{4}, [1, 3]))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeB{4}, [2, 4]))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeB{3}))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeB{4}))
    @test det_tangent_is_anticanonical(SGr(1, 4))
    @test det_tangent_is_anticanonical(LGr(2, 4))
    @test det_tangent_is_anticanonical(LGr(3, 6))
    @test det_tangent_is_anticanonical(LGr(4, 8))
    @test det_tangent_is_anticanonical(SGr(2, 6))
    @test det_tangent_is_anticanonical(SGr(1, 8))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeC{4}, [2, 4]))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeC{4}, [1, 3]))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeC{3}))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeC{4}))
    @test det_tangent_is_anticanonical(quadric(4))
    @test det_tangent_is_anticanonical(quadric(6))
    @test det_tangent_is_anticanonical(quadric(8))
    @test det_tangent_is_anticanonical(OGr(3, 6))
    @test det_tangent_is_anticanonical(OGr(4, 8))
    @test det_tangent_is_anticanonical(OGr(2, 8))
    @test det_tangent_is_anticanonical(OGr(3, 8))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeD{5}, [2, 4]))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeD{5}, [1, 5]))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeD{4}))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeD{5}))
    @test det_tangent_is_anticanonical(adjoint_variety(TypeG2))
    @test det_tangent_is_anticanonical(coadjoint_variety(TypeG2))
    @test det_tangent_is_anticanonical(full_flag_variety(TypeG2))
    @test det_tangent_is_anticanonical(adjoint_variety(TypeF4))
    @test det_tangent_is_anticanonical(coadjoint_variety(TypeF4))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeF4, [1, 3]))
    @test det_tangent_is_anticanonical(cayley_plane())
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeE{6}, 2))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeE{6}, [1, 3]))
    @test det_tangent_is_anticanonical(freudenthal_variety())
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeE{7}, 1))
    @test det_tangent_is_anticanonical(partial_flag_variety(TypeE{8}, 8))
  end

  @testset "Canonical and anticanonical bundles" begin
    # Rank is always 1
    @test rank_bundle(canonical_bundle(projective_space(3))) == 1
    @test rank_bundle(anticanonical_bundle(projective_space(3))) == 1

    # K_{P^1} = O(-2): H^1(P^1, O(-2)) = 1
    X1 = projective_space(1)
    K = canonical_bundle(X1)
    @test dimensions(K)[1] == 1
    @test dimensions(K)[0] == 0

    # -K_{P^1} = O(2): H^0(P^1, O(2)) = 3
    aK = anticanonical_bundle(X1)
    @test dimensions(aK)[0] == 3

    # K_{P^4} = O(-5): H^4(P^4, O(-5)) = 1, H^0 = 0
    X4 = projective_space(4)
    K4 = canonical_bundle(X4)
    @test dimensions(K4)[4] == 1
    @test dimensions(K4)[0] == 0

    # -K_{Gr(2,4)} = O(4): H^0(Gr(2,4), O(4)) = 105
    G24 = Gr(2, 4)
    aK24 = anticanonical_bundle(G24)
    @test rank_bundle(aK24) == 1
    @test dimensions(aK24)[0] == 105

    # canonical ⊗ anticanonical = structure sheaf (trivial bundle)
    X = projective_space(2)
    prod = canonical_bundle(X) ⊗ anticanonical_bundle(X)
    @test dimensions(prod)[0] == 1  # H^0(O) = 1
  end

  @testset "Fano index: ambient varieties" begin
    # P^n has Fano index n+1
    @test fano_index(projective_space(1)) == 2
    @test fano_index(projective_space(4)) == 5

    # Gr(k,n) has Fano index n
    @test fano_index(Gr(2, 5)) == 5
    @test fano_index(Gr(3, 6)) == 6
    @test fano_index(Gr(2, 4)) == 4

    # n-dimensional quadric has Fano index n
    @test fano_index(quadric(4)) == 4
    @test fano_index(quadric(3)) == 3

    # Exceptional varieties
    @test fano_index(cayley_plane()) == 12
    @test fano_index(freudenthal_variety()) == 18

    # Picard rank > 1: fano_index is gcd of anticanonical degrees
    @test fano_index(partial_flag_variety(TypeA{3}, (1, 3))) ==
      gcd(anticanonical_degrees(partial_flag_variety(TypeA{3}, (1, 3)))...)
  end

  @testset "Fano index: zero loci" begin
    # Z(O(k)) in P^n: fano_index = (n+1) - k
    X = projective_space(4)
    @test fano_index(zero_locus(line_bundle(X, 3))) == 2   # -K = O(2)
    @test fano_index(zero_locus(line_bundle(X, 5))) == 0   # CY: -K = O(0)

    # Z(O(1) ⊕ O(1)) in P^4: a codim-2 surface, fano_index = 5 - 2 = 3
    E2 = direct_sum(line_bundle(X, 1), line_bundle(X, 1))
    @test fano_index(zero_locus(E2)) == 3

    # Adjunction consistency: fano_index(Z) = fano_index(X) - deg(det E)
    X5 = Gr(2, 5)
    E5 = direct_sum(line_bundle(X5, 3), line_bundle(X5, 1))  # det = O(4)
    @test fano_index(zero_locus(E5)) == fano_index(X5) - 4

    # Throws for Picard rank > 1 ambient
    X_flag = partial_flag_variety(TypeA{3}, (1, 3))
    E_flag = structure_sheaf(X_flag)
    @test_throws ArgumentError fano_index(zero_locus(E_flag))
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

  @testset "Cohomology: higher-rank bundles" begin
    X = Gr(2, 5)
    S = universal_subbundle(X)
    Sym2 = symmetric_power(dual(S), 2)
    H = dimensions(cohomology(Sym2))
    @test H[0] == 15
    for i in 1:6
      @test H[i] == 0
    end
    @test euler_characteristic(Sym2) == 15

    X4 = Gr(2, 4)
    S4 = universal_subbundle(X4)
    Q4 = universal_quotient_bundle(X4)
    T = tensor_product(dual(S4), Q4)
    @test rank_bundle(T) == 4
    HT = dimensions(cohomology(T))
    @test HT[0] == 15
    for i in 1:4
      @test HT[i] == 0
    end

    X5 = Gr(2, 5)
    S5 = universal_subbundle(X5)
    W = exterior_power(dual(S5), 2)
    @test rank_bundle(W) == 1
    HW = dimensions(cohomology(W))
    @test HW[0] == 10
  end

  @testset "hilbert_polynomial" begin
    # χ(O_{P^3}(t)) = (t+1)(t+2)(t+3)/6 = 1 + 11t/6 + t^2 + t^3/6.
    P3 = projective_space(3)
    @test hilbert_polynomial(structure_sheaf(P3)) ==
      Rational{BigInt}[1, 11 // 6, 1, 1 // 6]

    # χ(O_{P^3}(1)(t)) = (t+2)(t+3)(t+4)/6 = 4 + 13t/3 + 3t^2/2 + t^3/6.
    @test hilbert_polynomial(line_bundle(P3, 1)) ==
      Rational{BigInt}[4, 13 // 3, 3 // 2, 1 // 6]

    # Polarization on Picard rank > 1: χ(O(t,t)) = (t+1)^2 on P^1 × P^1.
    X = projective_space(1) * projective_space(1)
    @test hilbert_polynomial(structure_sheaf(X), line_bundle(X, [1, 1])) ==
      Rational{BigInt}[1, 2, 1]

    # No polarization on Picard rank > 1: throws.
    @test_throws ArgumentError hilbert_polynomial(structure_sheaf(X))

    # Polarization must be a line bundle on the same variety as E.
    @test_throws ArgumentError hilbert_polynomial(structure_sheaf(P3), tangent_bundle(P3))
    @test_throws ArgumentError hilbert_polynomial(structure_sheaf(P3),
      line_bundle(projective_space(2), 1))
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
    V = flag_variety(4, [1, 2])
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

    V = LGr(3, 6)
    @test dimension(V) == 6
    @test LGr(3, 6) == SGr(3, 6)
    # The second argument is the ambient dimension and must equal 2n.
    @test_throws ArgumentError LGr(3, 7)
    @test_throws ArgumentError LGr(3, 5)
  end

  @testset "IGr (isotropic Grassmannian)" begin
    # IGr is a synonym of OGr, retained to emphasize the isotropy convention.
    @test IGr(2, 7) == OGr(2, 7)
    @test IGr(1, 9) == OGr(1, 9)
    @test IGr(3, 8) == OGr(3, 8)
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
    @test marked_dynkin_diagram(marked_dynkin_type(Gr(2, 5))) ==
      "○───×───○───○\n1   2   3   4"

    @test marked_dynkin_diagram(marked_dynkin_type(OGr(2, 8))) ==
      "        ○ 4\n       /\n○───×───○\n1   2   3"

    @test marked_dynkin_diagram(marked_dynkin_type(cayley_plane())) ==
      "        ○ 2\n        |\n×───○───○───○───○\n1   3   4   5   6"

    @test marked_dynkin_diagram(
      marked_dynkin_type(projective_space(2) * projective_space(2))
    ) == "A2:\n○───○\n1   2\n\nA2:\n○───○\n1   2\n(marked: 1, 3)"
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
      MDT = mdt(DT, marks)
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
    perm_A4 = levi_permutation(mdt(TypeA{4}, (2,)))
    @test collect(perm_A4) == [1, 2, 3]

    # D4/P1: unmarked nodes (2,3,4) of D4 form an A3 sub-diagram
    # cartan_type_with_ordering finds the canonical A3 ordering = [2,1,3]
    perm_D4 = levi_permutation(mdt(TypeD{4}, (1,)))
    @test collect(perm_D4) == [2, 1, 3]

    # B3/P1: unmarked nodes (2,3) form a B2 sub-diagram → identity
    perm_B3 = levi_permutation(mdt(TypeB{3}, (1,)))
    @test length(perm_B3) == 2

    # G2/P1: unmarked node (2) → trivial A1 sub-diagram
    perm_G2 = levi_permutation(mdt(TypeG2, (1,)))
    @test collect(perm_G2) == [1]
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  IrrepLevi round-trip for non-A Dynkin types
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "IrrepLevi round-trip: D types" begin
    MDT = mdt(TypeD{4}, (1,))
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
  #  Exterior power ranks
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
  #  H⁰ of exterior powers of tangent bundles
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "H⁰(∧ᵏT): A-types" begin
    T = tangent_bundle(Gr(2, 4))
    @test dimensions(exterior_power(T, 0))[0] == 1
    @test dimensions(exterior_power(T, 1))[0] == 15
    @test dimensions(exterior_power(T, 2))[0] == 90
    @test dimensions(exterior_power(T, 3))[0] == 175
    @test dimensions(exterior_power(T, 4))[0] == 105

    T2 = tangent_bundle(Gr(2, 5))
    @test dimensions(exterior_power(T2, 1))[0] == 24
    @test dimensions(exterior_power(T2, 2))[0] == 252
    @test dimensions(exterior_power(T2, 3))[0] == 1248
    @test dimensions(exterior_power(T2, 6))[0] == 1176

    T3 = tangent_bundle(Gr(3, 8))
    @test dimensions(exterior_power(T3, 1))[0] == 63
    @test dimensions(exterior_power(T3, 2))[0] == 1890
  end

  @testset "H⁰(∧ᵏT): B-types" begin
    T = tangent_bundle(partial_flag_variety(TypeB{3}, (1,)))
    @test dimensions(exterior_power(T, 1))[0] == 21
    @test dimensions(exterior_power(T, 2))[0] == 189
    @test dimensions(exterior_power(T, 3))[0] == 616
    @test dimensions(exterior_power(T, 4))[0] == 819
    @test dimensions(exterior_power(T, 5))[0] == 378
  end

  @testset "H⁰(∧ᵏT): C-types" begin
    T = tangent_bundle(partial_flag_variety(TypeC{3}, (1,)))
    @test dimensions(exterior_power(T, 1))[0] == 35
    @test dimensions(exterior_power(T, 2))[0] == 280
    @test dimensions(exterior_power(T, 3))[0] == 840
    @test dimensions(exterior_power(T, 4))[0] == 1050
    @test dimensions(exterior_power(T, 5))[0] == 462
  end

  @testset "H⁰(∧ᵏT): D-types" begin
    T = tangent_bundle(partial_flag_variety(TypeD{4}, (1,)))
    @test dimensions(exterior_power(T, 1))[0] == 28
    @test dimensions(exterior_power(T, 2))[0] == 350
    @test dimensions(exterior_power(T, 3))[0] == 1680
    @test dimensions(exterior_power(T, 6))[0] == 1386
  end

  @testset "H⁰(∧ᵏT): G2" begin
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
    @test variety(T) === X
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
    X1 = Gr(2, 5)
    @test rank_bundle(universal_subbundles(X1)[1]) == 2
    @test rank_bundle(universal_quotient_bundle(X1)) == 3

    # OGr(3, 7)
    X2 = OGr(3, 9)
    @test rank_bundle(universal_subbundle(X2)) == 3
    @test rank_bundle(residual_bundle(X2)) == 3

    # ℙ^4 = Gr(1, 5)
    X3 = projective_space(4)
    @test rank_bundle(universal_subbundle(X3)) == 1
    @test rank_bundle(universal_quotient_bundle(X3)) == 4

    # SGr(3,10)
    X4 = SGr(3, 10)
    U = universal_subbundle(X4)
    R = residual_bundle(X4)
    @test rank_bundle(U) == 3
    @test rank_bundle(R) == 4
    @test rank_bundle(universal_quotient_bundle(X4)) == 7
    @test dimensions(cohomology(dual(R) ⊗ U))[1] == 1

    #OGr(4,12)
    X5 = OGr(4, 12)
    U = universal_subbundle(X5, 1)
    R = residual_bundle(X5)
    @test rank_bundle(U) == 4
    @test rank_bundle(R) == 4
    @test rank_bundle(universal_quotient_bundle(X5)) == 8
    @test dimensions(cohomology(dual(R) ⊗ U))[1] == 1

    # Fl(1,2,4; 5)
    X6 = flag_variety(5, [1, 2, 4])
    τ = tautological_bundles(X6)
    Us = universal_subbundles(X6)
    @test length(τ) == 3
    @test rank_bundle(τ[1]) == 1
    @test rank_bundle(τ[2]) == 1
    @test rank_bundle(τ[3]) == 2
    @test dimensions(cohomology(dual(τ[3]) ⊗ τ[2]))[1] == 1
    @test dimensions(cohomology(dual(τ[2]) ⊗ τ[1]))[1] == 1
    @test rank_bundle(Us[1]) == 1
    @test rank_bundle(universal_subbundles(X6)[3]) == 4
    @test dimensions(cohomology(dual(τ[3]) ⊗ Us[2].pieces[2]))[1] == 1

    # universal_subbundle with explicit i > 1
    @test rank_bundle(universal_subbundle(X6, 2)) == 2
    @test rank_bundle(universal_subbundle(X6, 3)) == 4

    # Isotropic generalized Grassmannian: tautological_bundles = [U, R],
    # universal_subbundles = [U, U^⊥] (option (b)).
    @test rank_bundle.(tautological_bundles(OGr(2, 5))) == [2, 1]    # B_2/P_2
    @test rank_bundle.(universal_subbundles(OGr(2, 5))) == [2, 3]
    # The two-marked spinorial D_n/P_{n-1, n} also works (gives the natural
    # (n-1)-isotropic / unique containing max-isotropic flag).
    @test rank_bundle.(
      tautological_bundles(partial_flag_variety(TypeD{4}, (3, 4)))
    ) == [3, 1]
    @test rank_bundle.(
      universal_subbundles(partial_flag_variety(TypeD{4}, (3, 4)))
    ) == [3, 4]
    # Multi-step isotropic with D spinor-boundary at the last node.
    @test rank_bundle.(
      tautological_bundles(partial_flag_variety(TypeD{4}, (1, 3)))
    ) == [1, 3]    # m_2 = n-1, fix used: ω_n - ω_{n-1}
    # universal_subbundle(X, 2) on a B/C/D generalized Grassmannian returns U^⊥.
    @test rank_bundle(universal_subbundle(OGr(3, 9), 2)) == 6    # n-k = 9-3 = 6
    # Lagrangian / spinor cases: U^⊥ = U.
    @test rank_bundle(universal_subbundle(SGr(3, 6), 2)) ==
      rank_bundle(universal_subbundle(SGr(3, 6)))

    # residual_bundle: TypeB marked == R (OGr(3,7) = B_3/P_3, rank 7-2*3=1)
    @test rank_bundle(residual_bundle(OGr(3, 7))) == 1

    # residual_bundle: TypeB else, k < R-1 (OGr(2,9) = B_4/P_2, rank 9-2*2=5)
    @test rank_bundle(residual_bundle(OGr(2, 9))) == 5

    # residual_bundle: TypeC marked == R (SGr(4,8) = C_4/P_4, Lagrangian)
    @test rank_bundle(residual_bundle(SGr(4, 8))) == 0

    # residual_bundle: TypeD marked in (R, R-1) (spinor varieties)
    @test rank_bundle(residual_bundle(OGr(4, 10))) == 0  # D_5/P_4
    @test rank_bundle(residual_bundle(OGr(5, 10))) == 0  # D_5/P_5

    # residual_bundle: TypeD else, k < R-2 (OGr(2,10) = D_5/P_2, rank 10-2*2=6)
    @test rank_bundle(residual_bundle(OGr(2, 10))) == 6
  end

  @testset "Orthogonal & symplectic partial flag varieties" begin
    # tautological_bundles[i] is the i-th graded piece U_{m_i} / U_{m_{i-1}};
    # universal_subbundles[i] is the i-th subbundle U_{m_i}. In regular cases
    # rank(U_{m_i}) = m_i. The interesting case is type D at a spinor node
    # (m_i in {n-1, n}): U_{m_i} has rank n (max isotropic from one family).

    # --- Type B (OGr, odd quadric dim) ---
    # Generalized Grassmannians (1 marked node): [U, U^⊥].
    @test rank_bundle.(tautological_bundles(OGr(2, 7))) == [2, 3]    # B_3/P_2
    @test rank_bundle.(universal_subbundles(OGr(2, 7))) == [2, 5]
    @test rank_bundle.(tautological_bundles(OGr(3, 7))) == [3, 1]    # B_3/P_3 (max iso)
    @test rank_bundle.(universal_subbundles(OGr(3, 7))) == [3, 4]
    # Multi-step B partial flags.
    @test rank_bundle.(tautological_bundles(
      partial_flag_variety(TypeB{3}, (1, 2))
    )) == [1, 1]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeB{3}, (1, 2))
    )) == [1, 2]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeB{3}, (2, 3))     # last marked at R
    )) == [2, 3]
    # Full isotropic flag in B.
    @test rank_bundle.(
      universal_subbundles(
        partial_flag_variety(TypeB{4}, (1, 2, 3, 4))
      ),
    ) == [1, 2, 3, 4]

    # --- Type C (symplectic Grassmannian) ---
    @test rank_bundle.(tautological_bundles(SGr(2, 6))) == [2, 2]    # C_3/P_2
    @test rank_bundle.(universal_subbundles(SGr(2, 6))) == [2, 4]
    @test rank_bundle.(tautological_bundles(SGr(3, 6))) == [3, 0]    # Lagrangian: R = 0
    @test rank_bundle.(universal_subbundles(SGr(3, 6))) == [3, 3]    #   ⇒ U^⊥ = U
    # Multi-step C partial flags.
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeC{3}, (1, 2))
    )) == [1, 2]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeC{3}, (2, 3))     # last marked at R (Lagrangian step)
    )) == [2, 3]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeC{4}, (1, 2, 3))
    )) == [1, 2, 3]

    # --- Type D (even orthogonal) ---
    # Non-spinor generalized Grassmannian.
    @test rank_bundle.(tautological_bundles(OGr(2, 8))) == [2, 4]    # D_4/P_2, R rank 4
    @test rank_bundle.(universal_subbundles(OGr(2, 8))) == [2, 6]
    # Spinor varieties: U is max isotropic of dim n, R vanishes.
    @test rank_bundle.(tautological_bundles(OGr(3, 8))) == [4, 0]    # D_4/P_3 (- family)
    @test rank_bundle.(universal_subbundles(OGr(3, 8))) == [4, 4]
    @test rank_bundle.(tautological_bundles(OGr(4, 8))) == [4, 0]    # D_4/P_4 (+ family)
    # Multi-step D, no spinor node:
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{4}, (1, 2))
    )) == [1, 2]
    # Multi-step D, last marked = n-1 (the "minus" spinor): U_{n-1} has rank n.
    # The Bourbaki-Plate-IV correction ω_n - ω_{n-1} = L_n is required here;
    # the type-A formula ω_{m-1} - ω_m would give the wrong rank.
    @test rank_bundle.(tautological_bundles(
      partial_flag_variety(TypeD{4}, (1, 3))
    )) == [1, 3]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{4}, (1, 3))
    )) == [1, 4]
    # Multi-step D, last marked = n: U_n is max isotropic ("+" family).
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{4}, (1, 4))
    )) == [1, 4]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{4}, (2, 3))
    )) == [2, 4]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{4}, (2, 4))
    )) == [2, 4]
    # Two-marked spinorial D_n/P_{n-1, n} (= the (n-1)-isotropic Grassmannian,
    # Picard rank 2). U_{n-1} has rank n-1 here (genuine (n-1)-isotropic) and
    # U_n the unique containing max isotropic of rank n.
    @test rank_bundle.(tautological_bundles(
      partial_flag_variety(TypeD{4}, (3, 4))
    )) == [3, 1]
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{4}, (3, 4))
    )) == [3, 4]
    # 3-step D partial flag, last at n-1 (spinor).
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{5}, (1, 3, 4))
    )) == [1, 3, 5]
    # 3-step D partial flag, last at n.
    @test rank_bundle.(universal_subbundles(
      partial_flag_variety(TypeD{5}, (1, 3, 5))
    )) == [1, 3, 5]

    # --- universal_subbundle(X, i): index dispatch ---
    @test rank_bundle(universal_subbundle(OGr(2, 7), 2)) == 5         # B_3/P_2: U^⊥
    @test rank_bundle(universal_subbundle(OGr(2, 8), 2)) == 6         # D_4/P_2: U^⊥
    @test rank_bundle(universal_subbundle(
      partial_flag_variety(TypeD{4}, (1, 3)), 2
    )) == 4                                                            # spinor step
    @test_throws ArgumentError universal_subbundle(OGr(2, 7), 3)       # only 2 slots

    # --- Sanity: subs[k] always == sum of taut piece ranks (telescoping) ---
    for X in (
      OGr(2, 7), OGr(3, 7),
      partial_flag_variety(TypeB{3}, (1, 2)),
      partial_flag_variety(TypeB{4}, (1, 2, 3)),
      SGr(2, 6),
      partial_flag_variety(TypeC{3}, (1, 2)),
      partial_flag_variety(TypeD{4}, (1, 2)),
      partial_flag_variety(TypeD{4}, (1, 3)),
      partial_flag_variety(TypeD{4}, (3, 4)),
      partial_flag_variety(TypeD{5}, (1, 3, 4)),
    )
      taut = tautological_bundles(X)
      subs = universal_subbundles(X)
      @test rank_bundle(subs[end]) == sum(rank_bundle.(taut))
    end
  end

  @testset "is_exceptional_type" begin
    # Simple exceptional types
    @test is_exceptional_type(TypeE{6})
    @test is_exceptional_type(TypeE{7})
    @test is_exceptional_type(TypeE{8})
    @test is_exceptional_type(TypeF4)
    @test is_exceptional_type(TypeG2)

    # Simple classical types (negative cases)
    @test !is_exceptional_type(TypeA{4})
    @test !is_exceptional_type(TypeB{4})
    @test !is_exceptional_type(TypeC{4})
    @test !is_exceptional_type(TypeD{4})

    # Product types
    @test is_exceptional_type(ProductDynkinType{Tuple{TypeA{3},TypeG2}})
    @test !is_exceptional_type(ProductDynkinType{Tuple{TypeA{3},TypeB{4}}})

    # Variety level
    @test is_exceptional_type(partial_flag_variety("E6", 2))
    @test !is_exceptional_type(Gr(2, 5))
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

    # Q^3 = B_2/P_1: single spinor of rank 2
    X3 = quadric(3)
    S3 = spinor_bundle(X3)
    @test rank_bundle(S3) == 2

    # The two methods are mutually exclusive by type.
    # One-arg form is type B only: type D must pick a half-spinor instead.
    @test_throws ArgumentError spinor_bundle(X4)
    # Two-arg form is type D only: type B has a single spinor, use the one-arg form.
    @test_throws ArgumentError spinor_bundle(X3, :plus)
    @test_throws ArgumentError spinor_bundle(X3, :minus)
    # In type D the half symbol must be :plus or :minus.
    @test_throws ArgumentError spinor_bundle(X4, :other)
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

  @testset "print_hodge_diamond" begin
    # P^1: diamond is three rows (1 / 0 0 / 1); apex and base each carry a single "1".
    out1 = sprint(print_hodge_diamond, hodge_numbers(projective_space(1)))
    lines1 = filter(!isempty, split(out1, '\n'))
    @test length(lines1) == 3
    @test occursin("1", lines1[1])
    @test occursin("0", lines1[2])
    @test occursin("1", lines1[end])

    # Gr(2,4): 9 rows; h^{2,2} = b_4 = 2 sits on the middle row.
    out_gr = sprint(print_hodge_diamond, hodge_numbers(Gr(2, 4)))
    lines_gr = filter(!isempty, split(out_gr, '\n'))
    @test length(lines_gr) == 9
    @test occursin("2", lines_gr[5])
    @test occursin("1", lines_gr[1])
    @test occursin("1", lines_gr[end])
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

  @testset "HH^n via P[n]" begin
    P = hochschild_cohomology(projective_space(3))
    @test P[0] == 1
    @test P[1] == 15
    @test P[-1] == 0
    total = sum(P[n] for n in 0:6)
    expected = sum(P[p, q] for p in 0:3, q in 0:3)
    @test total == expected

    P_gr = hochschild_cohomology(Gr(2, 4))
    @test P_gr[0] == 1
    @test P_gr[1] == 15
  end

  @testset "PolyvectorParallelogram: show and Euler characteristic" begin
    P = hochschild_cohomology(projective_space(2))
    @test euler_characteristic(P) == 3

    P_gr = hochschild_cohomology(Gr(2, 4))
    @test euler_characteristic(P_gr) == 6

    out = sprint(show, MIME"text/plain"(), P)
    @test occursin("Polyvector parallelogram", out)
    @test occursin("dim = 2", out)

    @test sprint(show, P) == "PolyvectorParallelogram(dim=2)"
  end

  @testset "Twisted Hodge / Hochschild: exceptional types" begin
    let X = adjoint_variety(TypeG2)
      H = twisted_hodge_numbers(X, 0)
      @test H[1, 1] == 1
      @test H[1, 2] == 0
      @test H[2, 1] == 0
      d = dimension(X)
      for p in 0:d, q in 0:d
        p == q || @test H[p + 1, q + 1] == 0
      end

      P = hochschild_cohomology(X)
      @test P isa PolyvectorParallelogram
      @test P[0, 0] == 1
      @test P[1] == 14
    end

    let X = coadjoint_variety(TypeG2)
      H = twisted_hodge_numbers(X, 0)
      @test H[1, 1] == 1
      P = hochschild_cohomology(X)
      @test P[0, 0] == 1
    end
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

    var_counter = Ref(0)
    c_expr = PartialFlagVarieties.les_cokernel(
      [AffineExpr(1), AffineExpr(0)], [AffineExpr(2), AffineExpr(0)], var_counter
    )
    @test c_expr == [AffineExpr(1), AffineExpr(0)]
    @test var_counter[] == 0

    var_counter = Ref(0)
    c_zero = PartialFlagVarieties.les_cokernel(
      [AffineExpr(0), AffineExpr(0)], [AffineExpr(3), AffineExpr(5)], var_counter
    )
    @test c_zero == [AffineExpr(3), AffineExpr(5)]
    @test var_counter[] == 0

    # Case with connecting homomorphism ambiguity
    # H*(A) = [0,0,1], H*(B) = [0,0,0] on dim=2 variety
    a2 = Cohomology{BigInt}(BigInt[0, 0, 1], 2)
    b2 = Cohomology{BigInt}(BigInt[0, 0, 0], 2)
    (c2, det2) = solve_ses_cohomology(a2, b2)
    # χ(C) = χ(B) - χ(A) = 0 - 1 = -1
    @test sum((-1)^i * c2[i] for i in 0:2) == -1
  end

  @testset "Koszul symbolic solver" begin
    # AffineExpr arithmetic: a determined constant has no symbolic vars, and
    # adding/subtracting/scaling commutes with the constant/var split.
    e0 = AffineExpr(0)
    e3 = AffineExpr(3)
    @test is_determined(e0)
    @test is_determined(e3)
    @test is_zero_expr(e0)
    @test !is_zero_expr(e3)
    @test e3 + AffineExpr(2) == AffineExpr(5)
    @test e3 - e3 == AffineExpr(0)
    @test 2 * e3 == AffineExpr(6)

    # symbolic_variable(i) is x_i; subtracting it from itself eliminates the
    # variable, producing a fully determined zero expression.
    x1 = symbolic_variable(1)
    x2 = symbolic_variable(2)
    @test !is_determined(x1)
    @test !is_zero_expr(x1)
    @test is_zero_expr(x1 - x1)
    @test is_determined(x1 - x1)
    @test (x1 + x2) - x2 == x1
    @test 3 * x1 - x1 == 2 * x1

    # SES 0 → A → B → C → 0 with H*(A) = [1, 0], H*(B) = [2, 0] on a curve.
    # δ_0 ≤ a_1 = 0 forces δ_0 = 0, so C is fully determined and no symbolic
    # variable is introduced — the symbolic solver agrees with the numeric one.
    a = Cohomology{BigInt}(BigInt[1, 0], 1)
    b = Cohomology{BigInt}(BigInt[2, 0], 1)
    var_counter = Ref(0)
    c_sym = solve_ses_cohomology_symbolic(a, b, var_counter)
    @test c_sym isa Cohomology{AffineExpr}
    @test is_determined(c_sym[0])
    @test c_sym[0] == AffineExpr(1)
    @test c_sym[1] == AffineExpr(0)
    @test var_counter[] == 0

    # H*(A) = H*(B) = [1, 1]: the connecting map H^0(C) → H^1(A) has unknown
    # rank δ_0 ∈ {0, 1}, so the solver introduces a fresh symbolic variable.
    # The Euler characteristic is still pinned: χ(C) = χ(B) − χ(A) = 0.
    a_amb = Cohomology{BigInt}(BigInt[1, 1], 1)
    b_amb = Cohomology{BigInt}(BigInt[1, 1], 1)
    var_counter2 = Ref(0)
    c_amb = solve_ses_cohomology_symbolic(a_amb, b_amb, var_counter2)
    @test var_counter2[] >= 1
    @test c_amb[0] - c_amb[1] == AffineExpr(0)

    # Koszul filtration in P^3 for a line (CI of two hyperplanes), encoded as
    # the ambient cohomologies of the Koszul terms:
    #   K_0 = O, K_1 = O(−1)^2, K_2 = O(−2).
    # H*(P^3, K_i) = [0,0,0,0] for i ≥ 1 (negative twists in the BBW gap),
    # H*(P^3, K_0) = [1,0,0,0]. The Koszul filtration recovers
    # H*(line, O_line) = [1, 0], unambiguously, and no symbolic vars appear.
    koszul = [
      Cohomology{BigInt}(BigInt[1, 0, 0, 0], 3),
      Cohomology{BigInt}(BigInt[0, 0, 0, 0], 3),
      Cohomology{BigInt}(BigInt[0, 0, 0, 0], 3),
    ]
    var_counter3 = Ref(0)
    H_sym = solve_koszul_filtration_symbolic(koszul, 1, var_counter3)
    @test H_sym isa Cohomology{AffineExpr}
    @test H_sym[0] == AffineExpr(1)
    @test H_sym[1] == AffineExpr(0)
  end

  @testset "AffineExpr show" begin
    @test sprint(show, AffineExpr(5)) == "5"
    @test sprint(show, AffineExpr(0)) == "0"

    x1 = symbolic_variable(1)
    x2 = symbolic_variable(2)
    @test sprint(show, x1) == "x_1"
    @test sprint(show, -x1) == "-x_1"
    @test sprint(show, 2 * x1) == "2 * x_1"
    @test sprint(show, -3 * x1) == "-3 * x_1"
    @test sprint(show, AffineExpr(4) + x1) == "4 + x_1"
    @test sprint(show, x1 + x2) == "x_1 + x_2"
    @test sprint(show, x1 - x2) == "x_1 - x_2"
    @test sprint(show, 2 * x1 - x2) == "2 * x_1 - x_2"
    @test sprint(show, x1 - 2 * x2) == "x_1 - 2 * x_2"
  end

  @testset "long_exact_sequence_cokernel" begin
    # Numeric BigInt entry point: r=0 case wraps terms[1] in AffineExpr.
    var_counter = Ref(0)
    out0 = PartialFlagVarieties.long_exact_sequence_cokernel(
      Vector{BigInt}[BigInt[3, 5]], var_counter
    )
    @test out0 == AffineExpr[AffineExpr(3), AffineExpr(5)]
    @test var_counter[] == 0

    # r=1: two BigInt terms, delegates to les_cokernel.
    out1 = PartialFlagVarieties.long_exact_sequence_cokernel(
      Vector{BigInt}[BigInt[1, 0], BigInt[2, 0]], Ref(0)
    )
    @test out1 == AffineExpr[AffineExpr(1), AffineExpr(0)]

    # r=2: chains two LES solves.
    out2 = PartialFlagVarieties.long_exact_sequence_cokernel(
      Vector{BigInt}[BigInt[0, 0], BigInt[1, 0], BigInt[2, 0]], Ref(0)
    )
    @test length(out2) == 2

    # AffineExpr-valued entry point.
    affine_terms = [AffineExpr[AffineExpr(1), AffineExpr(0)],
      AffineExpr[AffineExpr(2), AffineExpr(0)]]
    out_aff = PartialFlagVarieties.long_exact_sequence_cokernel(affine_terms, Ref(0))
    @test out_aff == AffineExpr[AffineExpr(1), AffineExpr(0)]

    # r=0 in the AffineExpr method returns a copy of terms[1].
    single = [AffineExpr[AffineExpr(7), AffineExpr(0)]]
    out_single = PartialFlagVarieties.long_exact_sequence_cokernel(single, Ref(0))
    @test out_single == single[1]
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

  @testset "ZeroLocus: products" begin
    X12 = product(projective_space(1), projective_space(1))
    Z12 = zero_locus(line_bundle(X12, [1, 0]))
    Z3 = zero_locus(line_bundle(projective_space(1), 1))

    Z = product(Z12, Z3)
    ambient = product(projective_space(1), projective_space(1), projective_space(1))
    expected_bundle = direct_sum(
      line_bundle(ambient, [1, 0, 0]),
      line_bundle(ambient, [0, 0, 1]),
    )

    @test zerolocus62_label(Z) == zerolocus62_label(Z12 * Z3)
    @test ambient_variety(Z) == ambient
    @test dimension(Z) == dimension(Z12) + dimension(Z3)
    @test codimension(Z) == codimension(Z12) + codimension(Z3)
    @test rank_bundle(defining_bundle(Z)) ==
      rank_bundle(defining_bundle(Z12)) +
          rank_bundle(defining_bundle(Z3))
    @test zerolocus62_label(defining_bundle(Z)) == zerolocus62_label(expected_bundle)
  end

  # The zero locus of a section of O(1) on the Cayley plane OP² = E6/P1
  # is (up to isomorphism) the coadjoint variety of F4, i.e. F4/P4 = Coadj(F₄).
  # We verify this by checking that the Hilbert polynomial of Z (= h⁰(Z, O(t)))
  # matches that of F4/P4 for t = 0, 1, 2, 3 (which together uniquely determine the
  # polynomial of a 15-dimensional projective variety).
  @testset "ZeroLocus: O(1) on Cayley plane = coadjoint F4" begin
    X_E6 = cayley_plane()                         # E6/P1, dim = 16
    Y_F4 = coadjoint_variety(TypeF4)              # F4/P4, dim = 15
    Z = zero_locus(line_bundle(X_E6, 1))          # hyperplane section, dim = 15

    @test dimension(Z) == 15
    @test dimension(Z) == dimension(Y_F4)

    # Hilbert polynomial: h⁰(Z, O(t)) = dim H⁰(F4/P4, O(t)) for t = 0,1,2,3
    # t=0: both give 1 (structure sheaf)
    # t=1: both give 26 (sections of the hyperplane class)
    # t=2: both give 324
    # t=3: both give 2652
    for t in 0:3
      (H_Z, _) = cohomology_on_restriction(Z, line_bundle(X_E6, t))
      H_F4 = dimensions(cohomology(line_bundle(Y_F4, t)))
      @test H_Z[0] == H_F4[0]
    end

    # Test that cohomology_on_restriction validates that F lives on the ambient variety of Z
    F_wrong = line_bundle(Y_F4, 1)
    @test_throws ArgumentError cohomology_on_restriction(Z, F_wrong)
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
    @test is_calabi_yau(zero_locus(line_bundle(X, 5))) == true   # det = anticanonical
    @test is_calabi_yau(zero_locus(line_bundle(X, 3))) == false  # det ≠ anticanonical
    @test is_strict_calabi_yau(zero_locus(line_bundle(X, 5))) == true
    @test is_strict_calabi_yau(zero_locus(line_bundle(X, 3))) == false

    # Fano variety of lines is hyperkähler fourfold
    let X_gr = Gr(2, 6), S = universal_subbundle(X_gr)
      F = zero_locus(symmetric_power(dual(S), 3))
      @test is_calabi_yau(F) == true
      @test is_strict_calabi_yau(F) == false
    end

    X5 = projective_space(5)
    E = direct_sum(line_bundle(X5, 3), line_bundle(X5, 3))
    @test is_calabi_yau(zero_locus(E)) == true
    @test is_strict_calabi_yau(zero_locus(E)) == true
  end

  @testset "ZeroLocus: is_strongly_fano" begin
    P4 = projective_space(4)
    @test is_strongly_fano(zero_locus(line_bundle(P4, 3))) == true
    @test is_strongly_fano(zero_locus(line_bundle(P4, 4))) == true
    @test is_strongly_fano(zero_locus(line_bundle(P4, 5))) == false
    @test is_strongly_fano(zero_locus(line_bundle(P4, 6))) == false

    let X = Gr(2, 4), S = universal_subbundle(X)
      @test is_strongly_fano(zero_locus(symmetric_power(dual(S), 2))) == true
    end
    let X = Gr(2, 5), S = universal_subbundle(X)
      @test is_strongly_fano(zero_locus(symmetric_power(dual(S), 2))) == true
    end

    let X = projective_space(2) * projective_space(2)
      @test is_strongly_fano(zero_locus(line_bundle(X, [1, 1]))) == true
      @test is_strongly_fano(zero_locus(line_bundle(X, [3, 1]))) == false
    end
  end

  @testset "ZeroLocus: CY/Fano cross-type" begin
    Q5 = quadric(5)
    @test is_strongly_fano(zero_locus(line_bundle(Q5, 1))) == true
    @test is_strongly_fano(zero_locus(line_bundle(Q5, 5))) == false
    @test is_calabi_yau(zero_locus(line_bundle(Q5, 5))) == true

    S = SGr(2, 6)
    @test is_strongly_fano(zero_locus(line_bundle(S, 1))) == true

    D = OGr(2, 8)
    @test is_strongly_fano(zero_locus(line_bundle(D, 1))) == true
  end

  @testset "ZeroLocus: zero-dimensional" begin
    # Four generic hyperplanes in P^4 cut out a single reduced point.
    # A reduced point has trivial canonical (no dim > 0 entries to twist),
    # h^0(O) = 1, and -K is trivially ample: so CY, strict CY, and strongly Fano.
    P4 = projective_space(4)
    H = line_bundle(P4, 1)
    Z = zero_locus(reduce(direct_sum, [H for _ in 1:4]))
    @test dimension(Z) == 0
    @test euler_characteristic(Z) == 1
    @test is_calabi_yau(Z) == true
    @test is_strict_calabi_yau(Z) == true
    @test is_strongly_fano(Z) == true

    # A single hyperplane and one quadric in P^3 cut out two reduced points.
    # h^0(O_Z) = 2, so the locus is disconnected and not a strict CY.
    P3 = projective_space(3)
    Z2 = zero_locus(
      reduce(direct_sum, [line_bundle(P3, 1), line_bundle(P3, 1), line_bundle(P3, 2)])
    )
    @test dimension(Z2) == 0
    @test euler_characteristic(Z2) == 2
    @test is_calabi_yau(Z2) == false
    @test is_strict_calabi_yau(Z2) == false
  end

  @testset "ZeroLocus: koszul_terms" begin
    # Untwisted Koszul resolution of O_Z has rank(E) + 1 terms K_i = ∧^i E^∨.
    P4 = projective_space(4)
    Z = zero_locus(line_bundle(P4, 3))
    terms = koszul_terms(Z)
    @test length(terms) == 2
    @test terms[1] == structure_sheaf(P4)
    @test terms[2] == line_bundle(P4, -3)

    # rank(E) = 2 → 3 terms.
    E2 = direct_sum(line_bundle(P4, 1), line_bundle(P4, 2))
    terms2 = koszul_terms(zero_locus(E2))
    @test length(terms2) == 3

    # Twisted variant: each K_i is tensored by F.
    twist = line_bundle(P4, 1)
    terms_tw = koszul_terms(Z, twist)
    @test length(terms_tw) == 2
    @test terms_tw[1] == twist
    @test terms_tw[2] == tensor_product(twist, line_bundle(P4, -3))

    # F must live on the ambient.
    @test_throws ArgumentError koszul_terms(Z, structure_sheaf(Gr(2, 4)))
  end

  @testset "ZeroLocus: cohomology_on_restriction_symbolic" begin
    # Numeric-determined case (quintic CY3): same values as the numeric
    # restriction, wrapped in AffineExpr.
    P4 = projective_space(4)
    Z = zero_locus(line_bundle(P4, 5))
    H_sym = cohomology_on_restriction_symbolic(Z, Ref(0))
    @test H_sym isa Cohomology{AffineExpr}
    @test H_sym[0] == AffineExpr(1)
    @test H_sym[1] == AffineExpr(0)
    @test H_sym[2] == AffineExpr(0)
    @test H_sym[3] == AffineExpr(1)
    @test all(is_determined(H_sym[i]) for i in 0:3)

    # Twisted variant with explicit bundle argument.
    H_tw = cohomology_on_restriction_symbolic(Z, structure_sheaf(P4), Ref(0))
    @test H_tw[0] == AffineExpr(1)
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

  @testset "ZeroLocus: Mukai K3 surfaces" begin
    k3_hodge = [1 0 1; 0 20 0; 1 0 1]

    function test_k3_model(label, E)
      @testset "$label" begin
        Z = zero_locus(E)
        @test dimension(Z) == 2
        @test hodge_numbers(Z) == k3_hodge
      end
    end

    test_k3_model("g = 3: quartic in P^3", line_bundle(projective_space(3), 4))

    let X = projective_space(4)
      test_k3_model("g = 4: (2,3) in P^4", line_bundle(X, 2) + line_bundle(X, 3))
    end

    test_k3_model("g = 5: (2,2,2) in P^5", 3 * line_bundle(projective_space(5), 2))

    let X = Gr(2, 5)
      test_k3_model(
        "g = 6: (2,1,1,1) on Gr(2,5)", line_bundle(X, 2) + 3 * line_bundle(X, 1)
      )
    end

    let X = OGr(5, 10)
      # Mukai writes the generator as O(1/2); in this package it is line_bundle(X, 1).
      test_k3_model("g = 7: O(1/2)^8 on OGr+(5,10)", 8 * line_bundle(X, 1))
    end

    test_k3_model("g = 8: O(1)^6 on Gr(2,6)", 6 * line_bundle(Gr(2, 6), 1))

    let X = Gr(3, 6), U = universal_subbundle(X)
      test_k3_model("g = 9: ∧²U* + O(1)^4 on Gr(3,6)",
        exterior_power(dual(U), 2) + 4 * line_bundle(X, 1))
    end

    let X = Gr(2, 7), Q = universal_quotient_bundle(X)
      test_k3_model("g = 10: Q*(1) + O(1)^3 on Gr(2,7)",
        twist(dual(Q), 1) + 3 * line_bundle(X, 1))
    end

    let X = Gr(3, 7), U = universal_subbundle(X)
      test_k3_model("g = 12: (∧²U*)^3 + O(1) on Gr(3,7)",
        3 * exterior_power(dual(U), 2) + line_bundle(X, 1))
    end

    let X = Gr(3, 7), U = universal_subbundle(X), Q = universal_quotient_bundle(X)
      test_k3_model("g = 13: (∧²U*)^2 + ∧³Q on Gr(3,7)",
        2 * exterior_power(dual(U), 2) + exterior_power(Q, 3))
    end

    let X = OGr(3, 9)
      # This is the rank-2 irreducible bundle of highest weight ω₄.
      F = CompletelyReducibleBundle(X, [0, 0, 0, 1])
      test_k3_model("g = 18: F^5 on OGr(3,9)", 5 * F)
    end

    let X = Gr(4, 9), U = universal_subbundle(X)
      test_k3_model("g = 20: (∧²U*)^3 on Gr(4,9)", 3 * exterior_power(dual(U), 2))
    end
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

  @testset "Beilinson collection: type-A only" begin
    # The dual presentation Aₙ/Pₙ is still type A and must work.
    Xdual = Gr(3, 4)   # = dual ℙ³
    @test is_full_exceptional_sequence(beilinson_collection(Xdual), Xdual)
    @test is_full_exceptional_sequence(beilinson_collection_dual(Xdual), Xdual)

    # Beilinson is type-A only: the non-type-A models of ℙⁿ recognised by
    # is_projective_space must be rejected (their homogeneous bundles do not
    # reproduce the collection), rather than silently returning a wrong answer.
    for X in (SGr(1, 4), OGr(2, 5), OGr(2, 6), OGr(3, 6))  # C₂/P₁, B₂/P₂, D₃/P₂,₃
      @test is_projective_space(X)                # genuinely a ℙⁿ ...
      @test !(dynkin_type(X) <: TypeA)            # ... but not in type A
      @test_throws ArgumentError beilinson_collection(X)
      @test_throws ArgumentError beilinson_collection_dual(X)
    end
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

  @testset "Exceptional collection on zero locus: cubic threefold" begin
    # Cubic 3-fold Z(O(3)) ⊂ ℙ⁴ is Fano, O|_Z and O(1)|_Z are exceptional.
    X = projective_space(4)
    Z = zero_locus(line_bundle(X, 3))
    @test is_exceptional(structure_sheaf(X), Z)
    @test is_exceptional(line_bundle(X, 1), Z)
    @test is_exceptional_pair(structure_sheaf(X), line_bundle(X, 1), Z)
  end

  @testset "Exceptional sequence on zero locus: Fano lines on Q₁∩Q₂, g=3" begin
    # Conjecture: F₁(Q₁∩Q₂) for g=3 has (g-1)(2g-5) = 2 exceptional objects
    # F = Z((Sym²U^∨)^⊕2) on Gr(2,8)
    G = Gr(2, 8)
    U = universal_subbundle(G)
    E = 2 * symmetric_power(dual(U), 2)
    Z = zero_locus(E)
    @test dimension(Z) == dimension(G) - Int(rank_bundle(E))

    # Build the collection: for g=3, (g-1)(2g-5) = 2
    # k=0: Sym⁰U^∨ = O, Sym¹U^∨ = U^∨  (i=0,1 since g-1=2)
    # k≥g-3=0 means i < g-2=1, so for k=0 only i=0 survives in the short block
    # Actually for k < g-3: full block (g-1 objects); for k ≥ g-3: short block (g-2 objects)
    # g=3: g-3=0, so ALL k ≥ 0 use the short block of g-2=1 object
    # k ranges over 0..2g-6=0, so only k=0: just O_F
    # Wait, 2g-5=1 twist levels: k=0
    # Short block (k ≥ g-3=0): i = 0..g-3=0, so just O
    # That gives only 1 object — but (g-1)(2g-5)=2*1=2
    # Let me re-check: the conjecture says (g-3)(g-1) + (g-1)(g-2) objects
    # g=3: (0)(2) + (2)(1) = 2 objects
    # Full block for k=0..g-4=-1 (empty), short block for k=g-3..2g-5: k=0..1
    # Short block: i=0..g-3=0, so just {O(k)} for k=0,1
    L = CompletelyReducibleBundle[
      twist(structure_sheaf(G), 1, 0),  # O
      twist(structure_sheaf(G), 1, 1),  # O(1)
    ]
    @test length(L) == 2
    @test is_exceptional_sequence(L, Z)
    @test is_strong_exceptional_sequence(L, Z)
  end

  @testset "ZeroLocus: symbolic Hodge numbers" begin
    # Fully determined case: symbolic = numeric
    X = projective_space(4)
    Z = zero_locus(line_bundle(X, 5))
    H_sym = hodge_numbers_symbolic(Z)
    H_num = hodge_numbers(Z)
    for p in 0:3, q in 0:3
      @test is_determined(H_sym[p + 1, q + 1])
      @test is_determined(H_num[p + 1, q + 1])
      @test H_sym[p + 1, q + 1].constant == H_num[p + 1, q + 1].constant
    end

    # Constraint propagation: h^{1,0} = 0 for Küchle c6 Fano fourfold
    X2 = Gr(2, 6)
    E2 = reduce(direct_sum, [line_bundle(X2, 1) for _ in 1:4])
    Z2 = zero_locus(E2)
    H2 = hodge_numbers_symbolic(Z2)
    @test is_determined(H2[2, 1])       # h^{1,0} determined
    @test H2[2, 1].constant == 0        # h^{1,0} = 0
    @test is_determined(H2[2, 5])       # h^{1,4} determined
    @test H2[2, 5].constant == 0        # h^{1,4} = 0
    @test H2[2, 2].constant == 1        # h^{1,1} = 1
    @test H2[3, 3].constant == 8        # h^{2,2} = 8

    # Hodge + Serre symmetry in symbolic result
    X3 = projective_space(5)
    Z3 = zero_locus(line_bundle(X3, 3))
    H3 = hodge_numbers_symbolic(Z3)
    # h^{p,q} = h^{q,p} (Hodge symmetry)
    for p in 0:4, q in 0:4
      @test H3[p + 1, q + 1] == H3[q + 1, p + 1]
    end
    # h^{p,q} = h^{d-p,d-q} (Serre duality)
    for p in 0:4, q in 0:4
      @test H3[p + 1, q + 1] == H3[5 - p, 5 - q]
    end
  end

  @testset "ZeroLocus: Hodge numbers (HK fourfolds, K3^[2]-type)" begin
    # Fano variety of lines on a cubic fourfold — Beauville–Donagi (1985)
    # Zero locus of Sym³(S*) on Gr(2,6); dim = 8 - 4 = 4.
    # Hodge numbers of K3^[2]-type: h^{1,1}=21, h^{2,2}=232.
    X1 = Gr(2, 6)
    E1 = symmetric_power(dual(universal_subbundle(X1)), 3)
    Z1 = zero_locus(E1)
    @test dimension(Z1) == 4
    h1 = hodge_numbers(Z1)
    @test h1[1, 1] == 1   # h^{0,0}
    @test h1[2, 1] == 0   # h^{1,0}
    @test h1[2, 2] == 21  # h^{1,1}
    @test h1[2, 3] == 0   # h^{1,2}
    @test h1[2, 4] == 21  # h^{1,3}
    @test h1[2, 5] == 0   # h^{1,4}
    @test h1[3, 1] == 1   # h^{2,0}
    @test h1[3, 2] == 0   # h^{2,1}
    @test h1[3, 3] == 232 # h^{2,2}
    @test h1[3, 4] == 0   # h^{2,3}
    @test h1[3, 5] == 1   # h^{2,4}

    # Debarre–Voisin variety — Debarre–Voisin (2010)
    # Zero locus of ∧³(S*) on Gr(6,10); dim = 24 - 20 = 4.
    # Also of K3^[2]-type: same Hodge numbers.
    X2 = Gr(6, 10)
    E2 = exterior_power(dual(universal_subbundle(X2)), 3)
    Z2 = zero_locus(E2)
    @test dimension(Z2) == 4
    h2 = hodge_numbers(Z2)
    @test h2[1, 1] == 1   # h^{0,0}
    @test h2[2, 2] == 21  # h^{1,1}
    @test h2[3, 1] == 1   # h^{2,0}
    @test h2[3, 3] == 232 # h^{2,2}

    # Both Hodge diamonds agree
    @test h1 == h2
  end

  @testset "ZeroLocus: Hodge numbers (Küchle Fano fourfolds)" begin
    # Selected families from Küchle, Math. Z. 218 (1995), 563–575.
    # Invariants verified against the paper and by χ_top = 2 + 2b₂ + 2h¹³ + h²²
    # (using b₃=0 and b₂=h¹¹ for these Fano fourfolds).
    # Note: h^{p,q} is at matrix index h[p+1, q+1].
    # These families are Fano (not CY), so Serre duality H^k(F) = H^{d-k}(F*)
    # must NOT be applied when computing cohomology of bundles on the zero locus.

    # b2: O(2)² on Gr(2,5)   — h¹¹=1, h¹³=20, h²²=132, χ_top=176
    let X = Gr(2, 5), E = direct_sum(line_bundle(X, 2), line_bundle(X, 2)),
      Z = zero_locus(E), h = hodge_numbers(Z)

      @test dimension(Z) == 4
      @test h[2, 2] == 1     # h^{1,1}
      @test h[2, 4] == 20    # h^{1,3}
      @test h[3, 3] == 132   # h^{2,2}
      # χ_top = 2 + 2b₂ + 2h¹³ + h²²
      @test 2 + 2 * h[2, 2] + 2 * h[2, 4] + h[3, 3] == 176
    end

    # b6: O(1)³ + O(2) on Gr(2,6)   — h¹¹=1, h¹³=15, h²²=106, χ_top=140
    let X = Gr(2, 6), O1 = line_bundle(X, 1),
      E = direct_sum(direct_sum(O1, direct_sum(O1, O1)), line_bundle(X, 2)),
      Z = zero_locus(E), h = hodge_numbers(Z)

      @test dimension(Z) == 4
      @test h[2, 2] == 1     # h^{1,1}
      @test h[2, 4] == 15    # h^{1,3}
      @test h[3, 3] == 106   # h^{2,2}
      @test 2 + 2 * h[2, 2] + 2 * h[2, 4] + h[3, 3] == 140
    end

    # b7: O(1)⁶ on Gr(2,7)   — h¹¹=1, h¹³=6, h²²=57, χ_top=73
    let X = Gr(2, 7), O1 = line_bundle(X, 1),
      E = foldl(direct_sum, [line_bundle(X, 1) for _ in 1:6]), Z = zero_locus(E),
      h = hodge_numbers(Z)

      @test dimension(Z) == 4
      @test h[2, 2] == 1     # h^{1,1}
      @test h[2, 4] == 6     # h^{1,3}
      @test h[3, 3] == 57    # h^{2,2}
      @test 2 + 2 * h[2, 2] + 2 * h[2, 4] + h[3, 3] == 73
    end

    # b10: O(2) + Q*(1) on Gr(2,7)   — h¹¹=1, h¹³=14, h²²=100, χ_top=132
    # (Q*(1) = dual(universal_quotient_bundle) ⊗ O(1))
    let X = Gr(2, 7), Q = universal_quotient_bundle(X),
      E = direct_sum(line_bundle(X, 2), twist(dual(Q), 1)), Z = zero_locus(E),
      h = hodge_numbers(Z)

      @test dimension(Z) == 4
      @test h[2, 2] == 1     # h^{1,1}
      @test h[2, 4] == 14    # h^{1,3}
      @test h[3, 3] == 100   # h^{2,2}
      @test 2 + 2 * h[2, 2] + 2 * h[2, 4] + h[3, 3] == 132
    end

    # c3: Q*(1)² on Gr(3,7)   — h¹¹=1, h¹³=0, h²²=15, χ_top=19
    let X = Gr(3, 7), Q = universal_quotient_bundle(X), Qd1 = twist(dual(Q), 1),
      E = direct_sum(Qd1, Qd1), Z = zero_locus(E), h = hodge_numbers(Z)

      @test dimension(Z) == 4
      @test h[2, 2] == 1     # h^{1,1}
      @test h[2, 4] == 0     # h^{1,3}
      @test h[3, 3] == 15    # h^{2,2}
      @test 2 + 2 * h[2, 2] + 2 * h[2, 4] + h[3, 3] == 19
    end

    # c7: O(1) + (∧²Q* ⊗ O(1)) on Gr(3,8)   — h¹¹=2, h¹³=1, h²²=22, χ_top=30
    # This was previously catastrophically wrong (h²²=-1888) due to a spurious
    # Serre duality application on the Fano zero locus.
    let X = Gr(3, 8), Q = universal_quotient_bundle(X),
      E = direct_sum(line_bundle(X, 1), twist(exterior_power(dual(Q), 2), 1)),
      Z = zero_locus(E), h = hodge_numbers(Z)

      @test dimension(Z) == 4
      @test h[2, 2] == 2     # h^{1,1}
      @test h[2, 4] == 1     # h^{1,3}
      @test h[3, 3] == 22    # h^{2,2}
      @test 2 + 2 * h[2, 2] + 2 * h[2, 4] + h[3, 3] == 30
    end
  end

  # ════════════════════════════════════════════════════════════════════════════
  #  Zero locus: Hodge numbers dim=6 (2·Sym²U on Gr(2,8))
  # ════════════════════════════════════════════════════════════════════════════

  @testset "ZeroLocus: Hodge numbers dim=6" begin
    # 2·Sym²(U^∨) on Gr(2,8): a 6-dimensional zero locus.
    # Some entries are undetermined by Koszul + symmetry constraints;
    # those carry symbolic variables and are NOT defaulted to 0.
    let X = Gr(2, 8), U = universal_subbundle(X), E = 2 * symmetric_power(dual(U), 2),
      Z = zero_locus(E)

      @test dimension(Z) == 6
      @test euler_characteristic(Z) == 1

      h = hodge_numbers(Z)

      # Non-negative constant parts for all entries
      for p in 0:6, q in 0:6
        @test h[p + 1, q + 1].constant >= 0
      end

      # Determined entries
      @test is_determined(h[1, 1])  # h^{0,0}
      @test h[1, 1] == 1
      @test is_determined(h[7, 7])  # h^{6,6}
      @test h[7, 7] == 1
      @test is_determined(h[2, 2])  # h^{1,1}
      @test h[2, 2] == 1

      # h^{1,2} has constant part 3 (undetermined — symbolic expression)
      @test h[2, 3].constant == 3
      @test !is_determined(h[2, 3])

      # h^{2,5} = 0 (correctly enforced by combined Hodge–Serre constraint)
      @test is_determined(h[3, 6])
      @test h[3, 6] == 0

      # Symmetry: h^{p,q} = h^{d-p,d-q} (Serre duality)
      for p in 0:6, q in 0:6
        @test h[p + 1, q + 1] == h[7 - p, 7 - q]
      end
    end
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Zero loci: Hochschild cohomology (polyvector parallelogram)
  #  Reference values from Sage twisted-hodge-ci project (Brückmann formula)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "ZeroLocus: Hochschild cohomology" begin
    # Cubic surface (dim=2): degree 3 in ℙ³
    let X = projective_space(3), Z = zero_locus(line_bundle(X, 3))
      P = hochschild_cohomology(Z)
      @test P.dim == 2
      @test P[0, 0] == 1
      @test P[1, 1] == 4
      @test P[2, 0] == 4
      @test P[0, 1] == 0
      @test P[1, 0] == 0
      @test P[0, 2] == 0
    end

    # Quartic threefold (dim=3): Fano, ω = O(-1)
    let X = projective_space(4), Z = zero_locus(line_bundle(X, 4))
      P = hochschild_cohomology(Z)
      @test P.dim == 3
      @test P[0, 0] == 1
      @test P[1, 1] == 45
      @test P[2, 2] == 15
      @test P[3, 0] == 5
      # Vanishing checks
      @test P[1, 0] == 0
      @test P[2, 0] == 0
      @test P[2, 1] == 0
      @test P[3, 1] == 0
    end

    # Quintic CY3 (dim=3): Calabi–Yau, ω ≅ O
    let X = projective_space(4), Z = zero_locus(line_bundle(X, 5))
      P = hochschild_cohomology(Z)
      @test P.dim == 3
      @test P[0, 0] == 1
      @test P[0, 3] == 1
      @test P[1, 1] == 101
      @test P[1, 2] == 1
      @test P[2, 1] == 1
      @test P[2, 2] == 101
      @test P[3, 0] == 1
      @test P[3, 3] == 1
    end

    # Two quadrics in ℙ⁵ (dim=3)
    let X = projective_space(5), Z = zero_locus(2 * line_bundle(X, 2))
      P = hochschild_cohomology(Z)
      @test P.dim == 3
      @test P[0, 0] == 1
      @test P[1, 1].constant == 3  # undetermined connecting map rank
      @test P[2, 0] == 15
      @test P[3, 0] == 19
    end

    # CI(2,3) in ℙ⁴ — K3 surface
    let X = projective_space(4), Z = zero_locus(line_bundle(X, 2) + line_bundle(X, 3))
      P = hochschild_cohomology(Z)
      @test P.dim == 2
      @test P[0, 0] == 1
      @test P[0, 2] == 1
      @test P[1, 1] == 20
      @test P[2, 0] == 1
      @test P[2, 2] == 1
    end

    # CI(3,3) in ℙ⁵ — CY3
    let X = projective_space(5), Z = zero_locus(2 * line_bundle(X, 3))
      P = hochschild_cohomology(Z)
      @test P.dim == 3
      @test P[0, 0] == 1
      @test P[0, 3] == 1
      @test P[1, 1] == 73
      @test P[1, 2] == 1
      @test P[2, 1] == 1
      @test P[2, 2] == 73
      @test P[3, 0] == 1
      @test P[3, 3] == 1
    end
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  euler_characteristic_tangent_bundle  (issue #15)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "ZeroLocus: euler_characteristic_tangent_bundle" begin
    # ── Hypersurfaces in ℙ⁴ ─────────────────────────────────────────────────
    let X = projective_space(4)
      @test euler_characteristic_tangent_bundle(zero_locus(line_bundle(X, 3))) == -10
      @test euler_characteristic_tangent_bundle(zero_locus(line_bundle(X, 4))) == -45
      @test euler_characteristic_tangent_bundle(zero_locus(line_bundle(X, 5))) == -100  # CY3
    end

    # ── Hypersurfaces in ℙ⁵ ─────────────────────────────────────────────────
    let X = projective_space(5)
      @test euler_characteristic_tangent_bundle(zero_locus(line_bundle(X, 2))) == 15   # Q⁴
      @test euler_characteristic_tangent_bundle(zero_locus(line_bundle(X, 3))) == -20
      @test euler_characteristic_tangent_bundle(zero_locus(line_bundle(X, 4))) == -90
      @test euler_characteristic_tangent_bundle(zero_locus(line_bundle(X, 5))) == -216
    end

    # ── Küchle Fano 4-folds ──────────────────────────────────────────────────
    let X = Gr(2, 5), E = direct_sum(line_bundle(X, 2), line_bundle(X, 2))
      @test euler_characteristic_tangent_bundle(zero_locus(E)) == -72   # b2
    end

    let X = Gr(2, 7), E = foldl(direct_sum, [line_bundle(X, 1) for _ in 1:6])
      @test euler_characteristic_tangent_bundle(zero_locus(E)) == -42   # b7
    end

    let X = Gr(3, 7), Qd1 = twist(dual(universal_quotient_bundle(X)), 1)
      @test euler_characteristic_tangent_bundle(zero_locus(direct_sum(Qd1, Qd1))) == -18  # c3
    end
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  tangent_cohomology  (issue #15)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "les_kernel" begin
    vc = Ref(0)
    # H*(C) = 0 forces H*(A) = H*(B)
    @test PartialFlagVarieties.les_kernel(
      AffineExpr.([3, 1, 0]), AffineExpr.([0, 0, 0]), vc
    ) == AffineExpr.([3, 1, 0])

    # 0 → a₀ → 5 → 2 → a₁ → 0 → ⋯ gives a₀ = 3 + a₁, and a₂ = 0
    a = PartialFlagVarieties.les_kernel(
      AffineExpr.([5, 0, 0]), AffineExpr.([2, 0, 0]), vc
    )
    @test a[1] - a[2] == AffineExpr(3)
    @test is_zero_expr(a[3])
  end

  @testset "ZeroLocus: tangent_cohomology" begin
    # Quintic threefold: h¹(T) = 101 deformations, h²(T) = h¹(Ω¹) = 1.
    let Z = zero_locus(line_bundle(projective_space(4), 5))
      H = tangent_cohomology(Z)
      @test [H[i] for i in 0:3] == AffineExpr.([0, 101, 1, 0])
    end

    # K3 of degree 14 in Gr(2,6): h¹(T) = 20.
    let X = Gr(2, 6),
      Z = zero_locus(reduce(direct_sum, [line_bundle(X, 1) for _ in 1:6]))

      H = tangent_cohomology(Z)
      @test [H[i] for i in 0:2] == AffineExpr.([0, 20, 0])
    end

    # Cubic threefold: h⁰ - h¹ is pinned by χ even though the pair is open.
    let Z = zero_locus(line_bundle(projective_space(4), 3))
      H = tangent_cohomology(Z)
      alt = PartialFlagVarieties._alternating_sum(H.entries)
      @test alt == AffineExpr(euler_characteristic_tangent_bundle(Z))
      @test is_zero_expr(H[2]) && is_zero_expr(H[3])
    end
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  SpectralSequence of a filtered bundle
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "SpectralSequence" begin
    # Cominuscule ambient: one filtration step, trivially degenerate.
    let F = filtered_tangent_bundle(Gr(2, 4))
      S = spectral_sequence(F)
      @test does_E1_degenerate(S)
      H = cohomology(F)
      @test all(is_determined, H.entries)
      @test H[0] == 15  # H^0(T) = sl(4)
      @test all(is_zero_expr(H[q]) for q in 1:dimension(Gr(2, 4)))
    end

    X = SGr(2, 6)
    Omega = dual(filtered_tangent_bundle(X))
    @test n_filtration_steps(Omega) == 2

    # Ground truth: H^q(X, Ω^p) is the diagonal Betti table.  Whatever the
    # spectral sequence machinery declares determined must match it.
    betti = betti_numbers(X)
    for p in 0:3
      H = cohomology(exterior_power(Omega, p))
      for q in 0:dimension(X)
        if is_determined(H[q])
          @test H[q] == (p == q ? betti[p + 1] : 0)
        end
      end
    end

    # A genuinely non-degenerate case: Ω³(-1) has E₁ = (2, 1) in degrees
    # (4, 5) but true cohomology (1, 0); the result must stay symbolic and
    # retain the χ-level relation.
    W = tensor_product(exterior_power(Omega, 3), line_bundle(X, -1))
    S = spectral_sequence(W)
    @test !does_E1_degenerate(S)
    H = cohomology(W)
    @test !is_determined(H[4])
    @test H[4] - H[5] == AffineExpr(1)

    # The isotypical components partition the E₁ page.
    iso = isotypical_components(S)
    @test !isempty(iso)
    @test all(!isempty(E1_page(T)) for T in values(iso))
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Hodge numbers over ambients with filtered tangent bundle
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "hodge_numbers: filtered ambient SGr(2,6)" begin
    X = SGr(2, 6)
    betti = betti_numbers(X)

    # Hyperplane section (Fano 6-fold): Lefschetz forces h^{p,q} = h^{p,q}(X)
    # below the middle degree.
    h = hodge_numbers(zero_locus(line_bundle(X, 1)))
    for p in 0:6, q in 0:6
      p + q < 6 || continue
      @test h[p + 1, q + 1] == (p == q ? betti[p + 1] : 0)
    end
    @test h[4, 4] == 4

    # Fano threefold O(1)^4: h^{1,1} = 1 by Lefschetz and χ(Ω¹) = 4
    # forces h^{1,2} = 5.
    E = reduce(direct_sum, [line_bundle(X, 1) for _ in 1:4])
    h3 = hodge_numbers(zero_locus(E))
    @test h3[2, 2] == 1
    @test h3[2, 3] == 5
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Kodaira–Akizuki–Nakano vanishing  (issue #4)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "twisted_hodge_numbers: Kodaira-Akizuki-Nakano" begin
    Z = zero_locus(line_bundle(projective_space(4), 3))
    M = twisted_hodge_numbers(Z, 1)
    @test all(M[p + 1, q + 1] == 0 for p in 0:3, q in 0:3 if p + q > 3)
    M2 = twisted_hodge_numbers(Z, -1)
    @test all(M2[p + 1, q + 1] == 0 for p in 0:3, q in 0:3 if p + q < 3)
  end

  @testset "is_ample_line_bundle" begin
    X = Gr(2, 4)
    @test is_ample_line_bundle(line_bundle(X, 1))
    @test !is_ample_line_bundle(structure_sheaf(X))
    @test !is_ample_line_bundle(line_bundle(X, -1))
    @test !is_ample_line_bundle(universal_subbundle(X))
    X2 = partial_flag_variety(TypeA{3}, (1, 3))
    @test is_ample_line_bundle(line_bundle(X2, [1, 2]))
    @test !is_ample_line_bundle(line_bundle(X2, [1, 0]))
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  Lefschetz hyperplane injection  (issue #17)
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "hodge_numbers: Lefschetz injection" begin
    # Z(Q^*(1) ⊕ O(1)) ⊂ Gr(2,8) is an ample divisor in Z(Q^*(1)), so
    # h^{p,q} agrees with the latter below the middle degree; this resolves
    # the symbolic variable the long exact sequences leave open.
    X = Gr(2, 8)
    E = tensor_product(dual(universal_quotient_bundle(X)), line_bundle(X, 1))
    Z = zero_locus(direct_sum(E, line_bundle(X, 1)))
    h = hodge_numbers(Z)
    @test all(is_determined, h)
    @test h[3, 3] == 3  # h^{2,2} = h^{2,2}(Z(E)) = b_4(Gr(2,8)) by Lefschetz
    @test h[4, 4] == 3  # h^{3,3} via Serre duality
    @test h[3, 4] == 2  # middle-degree row
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  hodge_numbers vs hodge_numbers_symbolic agreement
  # ═══════════════════════════════════════════════════════════════════════════

  @testset "hodge_numbers_les == hodge_numbers_symbolic" begin
    let X = projective_space(4), Z = zero_locus(line_bundle(X, 5))
      @test hodge_numbers(Z) == hodge_numbers_symbolic(Z)
    end

    let X = projective_space(5), Z = zero_locus(2 * line_bundle(X, 2))
      @test hodge_numbers(Z) == hodge_numbers_symbolic(Z)
    end

    let X = Gr(2, 5), S = universal_subbundle(X), Z = zero_locus(symmetric_power(S, 2))
      @test hodge_numbers(Z) == hodge_numbers_symbolic(Z)
    end

    let X = Gr(2, 7), S = universal_subbundle(X), Z = zero_locus(2 * symmetric_power(S, 2))
      h1 = hodge_numbers(Z)
      h2 = hodge_numbers_symbolic(Z)
      d = dimension(Z)
      for p in 0:d, q in 0:d
        @test h1[p + 1, q + 1] == h2[p + 1, q + 1]
      end
    end

    let Z = zero_locus("44.70")
      expected = AffineExpr.([
        1 0 0 0 0
        0 2 0 0 0
        0 0 2 0 0
        0 0 0 2 0
        0 0 0 0 1
      ])
      @test hodge_numbers(Z) == expected
      @test hodge_numbers_symbolic(Z) == expected
      @test hodge_numbers_les(Z) == expected
    end

    let Z = zero_locus("2044.5m")
      x0 = symbolic_variable(0)
      expected = map(
        x -> x isa AffineExpr ? x : AffineExpr(x),
        [
          1 0 0 0 0
          0 3 + x0 x0 0 0
          0 x0 7 + 2 * x0 x0 0
          0 0 x0 3 + x0 0
          0 0 0 0 1
        ],
      )
      H = hodge_numbers_symbolic(Z)
      @test H == expected
      @test hodge_numbers_les(Z) == expected
      @test all(!is_determined(e) || e.constant >= 0 for e in H)
    end
  end

  # ═══════════════════════════════════════════════════════════════════════════
  #  ZeroLocus62 labels
  # ═══════════════════════════════════════════════════════════════════════════
  @testset "ZeroLocus62 labels" begin
    # ─── Ambient-only encoding against specification worked examples ──────
    @testset "ambient encoding" begin
      @test zerolocus62_label(projective_space(1)) == "1"
      @test zerolocus62_label(projective_space(3)) == "30"
      @test zerolocus62_label(Gr(2, 4)) == "31"
      @test zerolocus62_label(Gr(3, 6)) == "53"
      @test zerolocus62_label(partial_flag_variety(TypeA{3}, (1, 3))) == "34"

      # B-type
      @test zerolocus62_label(quadric(5)) == "H0"
      @test zerolocus62_label(full_flag_variety(TypeB{5})) == "JU"

      # D-type: OGr+(5,10) = D5/P5
      @test zerolocus62_label(partial_flag_variety(TypeD{5}, 5)) == "iF"

      # E-type: E7/P7 (Freudenthal variety)
      @test zerolocus62_label(freudenthal_variety()) == "u11"
    end

    # ─── Product varieties ────────────────────────────────────────────────
    @testset "product varieties" begin
      # P^1 × P^1
      let DT = ProductDynkinType{Tuple{TypeA{1},TypeA{1}}}
        X = partial_flag_variety(DT, (1, 2))
        @test zerolocus62_label(X) == "11"
      end

      # (P^1)^5
      let DT = ProductDynkinType{
          Tuple{
            ProductDynkinType{
              Tuple{
                ProductDynkinType{
                  Tuple{
                    ProductDynkinType{Tuple{TypeA{1},TypeA{1}}},
                    TypeA{1}},
                },
                TypeA{1}},
            },
            TypeA{1}},
        }
        X = partial_flag_variety(DT, Tuple(1:5))
        @test zerolocus62_label(X) == "11111"
      end
    end

    # ─── Bundle encoding against specification worked examples ────────────
    @testset "bundle encoding" begin
      # O(1) on P^1
      let X = projective_space(1)
        @test zerolocus62_label(zero_locus(line_bundle(X, 1))) == "1.0"
      end

      # O ⊕ O(1) on P^1 (rank exceeds dim, encode via bundle directly)
      let X = projective_space(1)
        E = direct_sum(structure_sheaf(X), line_bundle(X, 1))
        @test zerolocus62_label(E) == "1.x10"
      end

      # O(1) ⊕ O(1) on P^1 (rank exceeds dim, encode via bundle directly)
      let X = projective_space(1)
        E = direct_sum(line_bundle(X, 1), line_bundle(X, 1))
        @test zerolocus62_label(E) == "1.00"
      end

      # O(1) on P^3
      let X = projective_space(3)
        @test zerolocus62_label(zero_locus(line_bundle(X, 1))) == "30.0"
      end

      # O(1,0) ⊕ O(0,1) on P^1 × P^1
      let DT = ProductDynkinType{Tuple{TypeA{1},TypeA{1}}}
        X = partial_flag_variety(DT, (1, 2))
        L1 = line_bundle(X, [1, 0])
        L2 = line_bundle(X, [0, 1])
        E = direct_sum(L1, L2)
        @test zerolocus62_label(E) == "11.10"
      end

      # O(1,1) on P^1 × P^1
      let DT = ProductDynkinType{Tuple{TypeA{1},TypeA{1}}}
        X = partial_flag_variety(DT, (1, 2))
        L = line_bundle(X, [1, 1])
        @test zerolocus62_label(zero_locus(L)) == "11.E"
      end
    end

    # ─── Ambient decoding / round-trip ────────────────────────────────────
    @testset "ambient decoding" begin
      @test marked_dynkin_type(PartialFlagVariety("1")) ==
        marked_dynkin_type(projective_space(1))
      @test marked_dynkin_type(PartialFlagVariety("30")) ==
        marked_dynkin_type(projective_space(3))
      @test marked_dynkin_type(PartialFlagVariety("31")) ==
        marked_dynkin_type(Gr(2, 4))
      @test dimension(PartialFlagVariety("53")) == dimension(Gr(3, 6))
      @test dimension(PartialFlagVariety("H0")) == dimension(quadric(5))
    end

    # ─── Ambient round-trip identity ──────────────────────────────────────
    @testset "ambient round-trip" begin
      for X in [
        projective_space(1),
        projective_space(3),
        Gr(2, 4),
        Gr(3, 6),
        quadric(5),
        partial_flag_variety(TypeA{3}, (1, 3)),
        partial_flag_variety(TypeD{5}, 5),
        freudenthal_variety(),
      ]
        label = zerolocus62_label(X)
        X2 = PartialFlagVariety(label)
        @test zerolocus62_label(X2) == label
      end
    end

    # ─── Zero locus decoding ─────────────────────────────────────────────
    @testset "zero locus decoding" begin
      # O(1) on P^1
      let Z = zero_locus("1.0")
        @test dimension(ambient_variety(Z)) == 1
        @test rank_bundle(defining_bundle(Z)) == 1
      end

      # O(1) on P^3
      let Z = zero_locus("30.0")
        @test dimension(ambient_variety(Z)) == 3
        @test dimension(Z) == 2
      end

      # O(1,1) on P^1 × P^1
      let Z = zero_locus("11.E")
        @test dimension(ambient_variety(Z)) == 2
        @test rank_bundle(defining_bundle(Z)) == 1
      end
    end

    # ─── Zero locus round-trip ────────────────────────────────────────────
    @testset "zero locus round-trip" begin
      let X = projective_space(1)
        Z = zero_locus(line_bundle(X, 1))
        label = zerolocus62_label(Z)
        Z2 = zero_locus(label)
        @test zerolocus62_label(Z2) == label
      end

      let X = projective_space(3)
        Z = zero_locus(direct_sum(line_bundle(X, 1), line_bundle(X, 1)))
        label = zerolocus62_label(Z)
        Z2 = zero_locus(label)
        @test zerolocus62_label(Z2) == label
      end
    end

    # ─── Error handling ───────────────────────────────────────────────────
    @testset "error handling" begin
      @test_throws ArgumentError zero_locus("1")
      @test_throws ArgumentError zero_locus("")
    end
  end

  @testset "CacheConfig" begin
    info0 = PartialFlagVarieties.cache_info()
    @test info0 isa NamedTuple
    @test :tensor_product in keys(info0)
    @test :bwb_pair in keys(info0)
    @test :marked_dynkin in keys(info0)

    @test PartialFlagVarieties.clear_caches!() === nothing

    PartialFlagVarieties.configure_caches!(budget=64 * 1024 * 1024)
    info_small = PartialFlagVarieties.cache_info()
    @test info_small.tensor_product.maxsize > 0

    @test PartialFlagVarieties.configure_caches!() === nothing
  end
end  # @testset "PartialFlagVarieties.jl"
