# ═══════════════════════════════════════════════════════════════════════════════
#  Hyperkaehler.jl — Hodge numbers of hyperKähler fourfolds of K3^{[2]}-type
#
#  Two families of hyperKähler manifolds of K3^[2]-type are realised as zero
#  loci of regular sections of equivariant bundles on Grassmannians:
#
#  1. Fano variety of lines F(X)  — zero locus of  Sym³(S*)  on  Gr(2,6)
#     where S is the tautological rank-2 subbundle.
#     X = a cubic fourfold in ℙ⁵.  F(X) parametrises lines on X.
#     Dimension: dim Gr(2,6) - rank Sym³(S*) = 8 - 4 = 4.
#     Reference: Beauville–Donagi (1985).
#
#  2. Debarre–Voisin variety  — zero locus of  ∧³(S*)  on  Gr(6,10)
#     where S is the tautological rank-6 subbundle.
#     Dimension: dim Gr(6,10) - rank ∧³(S*) = 24 - 20 = 4.
#     Reference: Debarre–Voisin (2010).
#
#  Both are (irreducible) hyperKähler fourfolds of K3^[2]-type, whose Hodge
#  diamond is:
#
#              1
#            0   0
#          1   21  1
#        0   0   0   0
#      1   21  232  21  1
#        0   0   0   0
#          1   21  1
#            0   0
#              1
#
#  Usage:
#    julia --project=. examples/Hyperkaehler.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties

# ─── 1. Fano variety of lines F(X) on cubic fourfold ─────────────────────────

println("=" ^ 70)
println("1. Fano variety of lines  —  Gr(2,6),  E = Sym³(S*)")
println("=" ^ 70)
println()

X1 = Gr(2, 6)
E1 = symmetric_power(universal_subbundle(X1), 3)   # Sym³(S*), rank 4
Z1 = zero_locus(E1)
@assert dimension(Z1) == 4  "Expected dim Z = 4, got $(dimension(Z1))"

H1 = hodge_numbers(Z1)
println("Hodge diamond of F(cubic fourfold) in Gr(2,6) (d = 4):")
print_hodge_diamond(stdout, H1)
println()

@assert H1[1, 1] == 1   "h^{0,0} ≠ 1"
@assert H1[2, 2] == 21  "h^{1,1} ≠ 21"
@assert H1[3, 3] == 232 "h^{2,2} ≠ 232"

# ─── 2. Debarre–Voisin variety ────────────────────────────────────────────────

println("=" ^ 70)
println("2. Debarre–Voisin variety  —  Gr(6,10),  E = ∧³(S*)")
println("=" ^ 70)
println()

X2 = Gr(6, 10)
E2 = exterior_power(universal_subbundle(X2), 3)    # ∧³(S*), rank 20
Z2 = zero_locus(E2)
@assert dimension(Z2) == 4  "Expected dim Z = 4, got $(dimension(Z2))"

H2 = hodge_numbers(Z2)
println("Hodge diamond of Debarre–Voisin variety in Gr(6,10) (d = 4):")
print_hodge_diamond(stdout, H2)
println()

@assert H2[1, 1] == 1   "h^{0,0} ≠ 1"
@assert H2[2, 2] == 21  "h^{1,1} ≠ 21"
@assert H2[3, 3] == 232 "h^{2,2} ≠ 232"

# ─── Comparison ──────────────────────────────────────────────────────────────

println("Both Hodge diamonds match: ", H1 == H2)
println()
println("The Hodge numbers confirm that both are hyperKähler fourfolds of")
println("K3^[2]-type: h^{1,1} = 21, h^{2,2} = 232 = 22 + 2*binom(22,2)/22.")
