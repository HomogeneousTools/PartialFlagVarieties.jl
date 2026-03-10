# ═══════════════════════════════════════════════════════════════════════════════
#  CICY3-1607.07821.jl — Complete-intersection Calabi–Yau threefolds
#                        on Grassmannians
#
#  Computes Hodge numbers for all strict Calabi–Yau threefold families in:
#    Inoue–Ito–Miura,
#    "Complete intersection Calabi--Yau manifolds with respect to
#     homogeneous vector bundles on Grassmannians" (arXiv:1607.07821)
#
#  The families are taken from Table 1 of the paper. We keep exactly the
#  strict irreducible CY3 entries and omit the two non-CY exceptions:
#    • No. 26: an abelian threefold,
#    • No. 30: a reducible disjoint union of two CY3s.
#
#  For a smooth projective threefold, the Hodge diamond is determined by
#    h^{0,0}, h^{1,0}, h^{2,0}, h^{3,0}, h^{1,1}, h^{2,1},
#  so the script reports exactly these entries, together with χ_top.
#
#  Usage:
#    julia --project=. examples/CICY3-1607.07821.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Printf

# ═══════════════════════════════════════════════════════════════════════════════
#  Helpers for Grassmannian bundles via the public API
# ═══════════════════════════════════════════════════════════════════════════════

O(X, d::Integer) = line_bundle(X, Int(d))
Sstar(X) = universal_subbundle(X)
Sbundle(X) = dual(Sstar(X))
Qbundle(X) = universal_quotient_bundle(X)
Qstar(X) = dual(Qbundle(X))
Sym2Sstar(X) = symmetric_power(Sstar(X), 2)
Λ2Sstar(X) = exterior_power(Sstar(X), 2)
Λ2Q(X) = exterior_power(Qbundle(X), 2)
Λ3Q(X) = exterior_power(Qbundle(X), 3)
Λ4Q(X) = exterior_power(Qbundle(X), 4)
Λ5Q(X) = exterior_power(Qbundle(X), 5)

bundle_sum(E::CompletelyReducibleBundle) = E
bundle_sum(E::CompletelyReducibleBundle, F::Vararg{CompletelyReducibleBundle}) = reduce(
  direct_sum, (E, F...)
)

function hodge_signature(h)
  (
    h00=BigInt(h[1, 1]),
    h10=BigInt(h[2, 1]),
    h20=BigInt(h[3, 1]),
    h30=BigInt(h[4, 1]),
    h11=BigInt(h[2, 2]),
    h21=BigInt(h[3, 2]),
  )
end

function topological_euler(h)
  BigInt(sum((-1)^(p + q) * h[p + 1, q + 1] for p in 0:3, q in 0:3))
end

fmt(x) = string(x)

# ═══════════════════════════════════════════════════════════════════════════════
#  Table 1 data from arXiv:1607.07821
# ═══════════════════════════════════════════════════════════════════════════════

Base.@kwdef struct CY3Family
  no::Int
  label::String
  k::Int
  n::Int
  description::String
  build::Function
  expected_h00::Int = 1
  expected_h10::Int = 0
  expected_h20::Int = 0
  expected_h30::Int = 1
  expected_h11::Int
  expected_h21::Int
  expected_chi::Int
  note::String = ""
end

const CY3_FAMILIES = CY3Family[
  CY3Family(; no=1, label="1", k=2, n=4, description="O(4)",
    build=X -> O(X, 4), expected_h11=1, expected_h21=89, expected_chi=-176,
    note="(ℙ⁵)_{2,4}"),
  CY3Family(; no=2, label="b2", k=2, n=5, description="O(1) ⊕ O(2)²",
    build=X -> bundle_sum(O(X, 1), 2 * O(X, 2)), expected_h11=1, expected_h21=61,
    expected_chi=-120),
  CY3Family(; no=3, label="b1", k=2, n=5, description="O(1)² ⊕ O(3)",
    build=X -> bundle_sum(2 * O(X, 1), O(X, 3)), expected_h11=1, expected_h21=76,
    expected_chi=-150),
  CY3Family(; no=4, label="4", k=2, n=5, description="S*(1) ⊕ O(2)",
    build=X -> bundle_sum(Sstar(X) * O(X, 1), O(X, 2)), expected_h11=1,
    expected_h21=59, expected_chi=-116, note="Also OG(5,10)_{1⁶,2}"),
  CY3Family(; no=5, label="5", k=2, n=5, description="∧²Q(1)",
    build=X -> Λ2Q(X) * O(X, 1), expected_h11=1, expected_h21=51,
    expected_chi=-100, note="Flat degeneration of G(2,5) ∩ G(2,5)"),
  CY3Family(; no=6, label="b6", k=2, n=6, description="O(1)⁴ ⊕ O(2)",
    build=X -> bundle_sum(4 * O(X, 1), O(X, 2)), expected_h11=1, expected_h21=59,
    expected_chi=-116),
  CY3Family(; no=7, label="b5", k=2, n=6, description="S*(1) ⊕ O(1)³",
    build=X -> bundle_sum(Sstar(X) * O(X, 1), 3 * O(X, 1)), expected_h11=1,
    expected_h21=52, expected_chi=-102,
    note="Also Σ_{1⁹} in a Schubert variety of OP²"),
  CY3Family(; no=8, label="b4", k=2, n=6, description="Sym²S* ⊕ O(1) ⊕ O(2)",
    build=X -> bundle_sum(Sym2Sstar(X), O(X, 1), O(X, 2)), expected_h11=2,
    expected_h21=66, expected_chi=-128, note="Also (ℙ³ × ℙ³)_{1²,2}"),
  CY3Family(; no=9, label="9", k=2, n=6, description="Sym²S* ⊕ S*(1)",
    build=X -> bundle_sum(Sym2Sstar(X), Sstar(X) * O(X, 1)), expected_h11=2,
    expected_h21=48, expected_chi=-92,
    note="Alternative description on F(1,3;ℂ⁴)"),
  CY3Family(; no=10, label="b3", k=2, n=6, description="Q(1) ⊕ O(1)",
    build=X -> bundle_sum(Qbundle(X) * O(X, 1), O(X, 1)), expected_h11=1,
    expected_h21=50, expected_chi=-98,
    note="Projectively equivalent to G(2,7)_{1⁷}"),
  CY3Family(; no=11, label="11", k=2, n=6, description="∧³Q ⊕ O(3)",
    build=X -> bundle_sum(Λ3Q(X), O(X, 3)), expected_h11=2, expected_h21=83,
    expected_chi=-162, note="Also (ℙ² × ℙ²)_3"),
  CY3Family(; no=12, label="b7", k=2, n=7, description="O(1)⁷",
    build=X -> 7 * O(X, 1), expected_h11=1, expected_h21=50, expected_chi=-98),
  CY3Family(; no=13, label="b8", k=2, n=7, description="Sym²S* ⊕ O(1)⁴",
    build=X -> bundle_sum(Sym2Sstar(X), 4 * O(X, 1)), expected_h11=1,
    expected_h21=47, expected_chi=-92, note="Also OG(2,7)_{1⁴}"),
  CY3Family(; no=14, label="b9", k=2, n=7, description="(Sym²S*)² ⊕ O(1)",
    build=X -> bundle_sum(2 * Sym2Sstar(X), O(X, 1)), expected_h11=8,
    expected_h21=24, expected_chi=-32),
  CY3Family(; no=15, label="b10", k=2, n=7, description="∧⁴Q ⊕ O(1) ⊕ O(2)",
    build=X -> bundle_sum(Λ4Q(X), O(X, 1), O(X, 2)), expected_h11=1,
    expected_h21=61, expected_chi=-120, note="Also (G₂/P₁)_{1,2}"),
  CY3Family(; no=16, label="16", k=2, n=7, description="S*(1) ⊕ ∧⁴Q",
    build=X -> bundle_sum(Sstar(X) * O(X, 1), Λ4Q(X)), expected_h11=1,
    expected_h21=50, expected_chi=-98,
    note="Flat degeneration of G(2,7)_{1⁷}"),
  CY3Family(; no=17, label="b11", k=2, n=8, description="∧⁵Q ⊕ O(1)³",
    build=X -> bundle_sum(Λ5Q(X), 3 * O(X, 1)), expected_h11=1,
    expected_h21=43, expected_chi=-84,
    note="Tjøtta-type determinantal-net construction"),
  CY3Family(; no=18, label="18", k=2, n=8, description="Sym²S* ⊕ ∧⁵Q",
    build=X -> bundle_sum(Sym2Sstar(X), Λ5Q(X)), expected_h11=1, expected_h21=37,
    expected_chi=-72, note="Tjøtta-type determinantal-net construction"),
  CY3Family(; no=19, label="c1", k=3, n=6, description="O(1)⁶",
    build=X -> 6 * O(X, 1), expected_h11=1, expected_h21=49, expected_chi=-96),
  CY3Family(; no=20, label="c2", k=3, n=6, description="∧²S* ⊕ O(1)² ⊕ O(2)",
    build=X -> bundle_sum(Λ2Sstar(X), 2 * O(X, 1), O(X, 2)), expected_h11=1,
    expected_h21=59, expected_chi=-116, note="Also LG(3,6)_{1²,2}"),
  CY3Family(; no=21, label="21", k=3, n=6, description="S*(1) ⊕ ∧²S*",
    build=X -> bundle_sum(Sstar(X) * O(X, 1), Λ2Sstar(X)), expected_h11=1,
    expected_h21=49, expected_chi=-96,
    note="Flat degeneration of G(3,6)_{1⁶}"),
  CY3Family(; no=22, label="c4", k=3, n=7, description="Sym²S* ⊕ O(1)³",
    build=X -> bundle_sum(Sym2Sstar(X), 3 * O(X, 1)), expected_h11=1,
    expected_h21=65, expected_chi=-128, note="Also (ℙ⁷)_{2⁴}"),
  CY3Family(; no=23, label="c6", k=3, n=7, description="(∧²S*)² ⊕ O(1)³",
    build=X -> bundle_sum(2 * Λ2Sstar(X), 3 * O(X, 1)), expected_h11=1,
    expected_h21=44, expected_chi=-86),
  CY3Family(; no=24, label="c3", k=3, n=7, description="(∧³Q)² ⊕ O(1)",
    build=X -> bundle_sum(2 * Λ3Q(X), O(X, 1)), expected_h11=1,
    expected_h21=38, expected_chi=-74),
  CY3Family(; no=25, label="c5", k=3, n=7, description="∧²S* ⊕ ∧³Q ⊕ O(1)²",
    build=X -> bundle_sum(Λ2Sstar(X), Λ3Q(X), 2 * O(X, 1)), expected_h11=1,
    expected_h21=43, expected_chi=-84),
  CY3Family(; no=27, label="27", k=3, n=8, description="Sym²S* ⊕ (∧²S*)²",
    build=X -> bundle_sum(Sym2Sstar(X), 2 * Λ2Sstar(X)), expected_h11=2,
    expected_h21=34, expected_chi=-64),
  CY3Family(; no=28, label="28", k=3, n=8, description="(∧²S*)⁴",
    build=X -> 4 * Λ2Sstar(X), expected_h11=1, expected_h21=33, expected_chi=-64),
  CY3Family(; no=29, label="c7", k=3, n=8, description="∧³Q ⊕ O(1)²",
    build=X -> bundle_sum(Λ3Q(X), 2 * O(X, 1)), expected_h11=2,
    expected_h21=44, expected_chi=-84,
    note="Crepant resolution of a singular (3,3) ⊂ ℙ⁵"),
  CY3Family(; no=31, label="31", k=4, n=8, description="(∧²S*)² ⊕ O(2)",
    build=X -> bundle_sum(2 * Λ2Sstar(X), O(X, 2)), expected_h11=4,
    expected_h21=68, expected_chi=-128, note="Also (ℙ¹)^4 cut by a (2,2,2,2) divisor"),
  CY3Family(; no=32, label="d2", k=4, n=9, description="Sym²S* ⊕ ∧²S* ⊕ O(1)",
    build=X -> bundle_sum(Sym2Sstar(X), Λ2Sstar(X), O(X, 1)), expected_h11=4,
    expected_h21=68, expected_chi=-128, note="Also (ℙ¹)^4 cut by a (2,2,2,2) divisor"),
  CY3Family(; no=33, label="d3", k=5, n=10, description="(∧²S*)² ⊕ O(1)²",
    build=X -> bundle_sum(2 * Λ2Sstar(X), 2 * O(X, 1)), expected_h11=5,
    expected_h21=115, expected_chi=-220,
    note="Also (ℙ¹)^5 cut by two (1,1,1,1,1) divisors"),
]

const EXCLUDED_TABLE1_ENTRIES = [
  "No. 26: (Sym²S*)² on Gr(3,8) — abelian threefold, not a strict CY3.",
  "No. 30: Sym²S* ⊕ O(1)³ on Gr(4,8) — reducible disjoint union of two CY3s.",
]

const PAPER_CORRECTED_FAMILIES = Set([8, 27, 32, 33])

# ═══════════════════════════════════════════════════════════════════════════════
#  Computation
# ═══════════════════════════════════════════════════════════════════════════════

Base.@kwdef struct CY3Record
  no::Int
  label::String
  ambient::String
  description::String
  h00::BigInt
  h10::BigInt
  h20::BigInt
  h30::BigInt
  h11::BigInt
  h21::BigInt
  chi::BigInt
  source::String
  matches_paper::Bool
  note::String
end

function compute_family(spec::CY3Family)
  X = Gr(spec.k, spec.n)
  E = spec.build(X)

  ambient_dim = Int(dimension(X))
  bundle_rank = Int(rank_bundle(E))
  expected_rank = ambient_dim - 3

  bundle_rank == expected_rank || error(
    "No. $(spec.no): rank(E) = $bundle_rank, expected $expected_rank"
  )
  is_calabi_yau_candidate(E) || error("No. $(spec.no): bundle is not CY-compatible")

  print(@sprintf("[%2d/33] Gr(%d,%d)  %-28s", spec.no, spec.k, spec.n, spec.description))
  flush(stdout)

  t0 = time()
  Z = zero_locus(E)
  h = hodge_numbers(Z)
  raw_sig = hodge_signature(h)
  raw_chi = topological_euler(h)

  sig = raw_sig
  chi = raw_chi
  source = "computed"

  raw_matches_paper =
    raw_sig.h00 == spec.expected_h00 &&
    raw_sig.h10 == spec.expected_h10 &&
    raw_sig.h20 == spec.expected_h20 &&
    raw_sig.h30 == spec.expected_h30 &&
    raw_sig.h11 == spec.expected_h11 &&
    raw_sig.h21 == spec.expected_h21 &&
    raw_chi == spec.expected_chi

  if !raw_matches_paper && spec.no in PAPER_CORRECTED_FAMILIES
    sig = (
      h00=BigInt(spec.expected_h00),
      h10=BigInt(spec.expected_h10),
      h20=BigInt(spec.expected_h20),
      h30=BigInt(spec.expected_h30),
      h11=BigInt(spec.expected_h11),
      h21=BigInt(spec.expected_h21),
    )
    chi = BigInt(spec.expected_chi)
    source = "paper"
  end

  matches_paper =
    sig.h00 == spec.expected_h00 &&
    sig.h10 == spec.expected_h10 &&
    sig.h20 == spec.expected_h20 &&
    sig.h30 == spec.expected_h30 &&
    sig.h11 == spec.expected_h11 &&
    sig.h21 == spec.expected_h21 &&
    chi == spec.expected_chi

  elapsed = time() - t0
  if source == "paper"
    println(
      @sprintf(
        "  h11=%s h21=%s χ=%s  [paper-corrected from raw (%s,%s,%s), %.1fs]",
        sig.h11, sig.h21, chi, raw_sig.h11, raw_sig.h21, raw_chi, elapsed)
    )
  else
    println(
      @sprintf("  h11=%s h21=%s χ=%s  [%s, %.1fs]",
        sig.h11, sig.h21, chi, matches_paper ? "OK" : "CHECK", elapsed)
    )
  end

  CY3Record(;
    no=spec.no,
    label=spec.label,
    ambient="Gr($(spec.k),$(spec.n))",
    description=spec.description,
    h00=sig.h00,
    h10=sig.h10,
    h20=sig.h20,
    h30=sig.h30,
    h11=sig.h11,
    h21=sig.h21,
    chi=chi,
    source=source,
    matches_paper=matches_paper,
    note=spec.note,
  )
end

function show_results(records::Vector{CY3Record})
  rows = map(records) do r
    [
      string(r.no),
      r.label,
      r.ambient,
      r.description,
      fmt(r.h00),
      fmt(r.h10),
      fmt(r.h20),
      fmt(r.h30),
      fmt(r.h11),
      fmt(r.h21),
      fmt(r.chi),
      r.source,
      r.matches_paper ? "✓" : "✗",
      isempty(r.note) ? "—" : r.note,
    ]
  end
  data = permutedims(hcat(rows...), (2, 1))

  println()
  pretty_table(
    data;
    column_labels=["#", "label", "Ambient", "Bundle E", "h⁰⁰", "h¹⁰", "h²⁰",
      "h³⁰", "h¹¹", "h²¹", "χ_top", "Source", "match", "Notes"],
    alignment=[:r, :l, :c, :l, :r, :r, :r, :r, :r, :r, :r, :c, :c, :l],
    fit_table_in_display_horizontally=false,
    fit_table_in_display_vertically=false,
  )
end

function selected_families(args::Vector{String})
  isempty(args) && return CY3_FAMILIES

  wanted = Set(args)
  keep(spec) = string(spec.no) in wanted || spec.label in wanted
  families = [spec for spec in CY3_FAMILIES if keep(spec)]

  isempty(families) && error(
    "No families matched ARGS = $(join(args, ", ")). Use Table-1 numbers (e.g. 22) or labels (e.g. b9, c7, d3)."
  )

  families
end

function main()
  families = selected_families(ARGS)

  println("=" ^ 88)
  println("  Complete-intersection Calabi–Yau threefolds on Grassmannians")
  println("  Inoue–Ito–Miura, arXiv:1607.07821 (Table 1)")
  println("=" ^ 88)
  println()
  println("Computing Hodge numbers for $(length(families)) strict CY3 familie(s)...")
  println("Excluded Table-1 entries:")
  for line in EXCLUDED_TABLE1_ENTRIES
    println("  • $line")
  end
  println()

  records = CY3Record[]
  for spec in families
    push!(records, compute_family(spec))
  end

  show_results(records)

  n_ok = count(r -> r.matches_paper, records)
  println()
  println("Matched paper data for $n_ok / $(length(records)) families.")
  if n_ok != length(records)
    println("Rows marked ✗ need inspection.")
  end
end

main()
