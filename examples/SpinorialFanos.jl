# ═══════════════════════════════════════════════════════════════════════════════
#  SpinorialFanos.jl — Spinorial Fano manifolds (Frassineti–Manivel)
#
#  Reproduces Hodge data for Table 2 of
#    A. Frassineti, L. Manivel, "Spinorial Fano manifolds", arXiv:2605.28712.
#
#  Each variety X is the smooth zero locus of a general section of a spinor
#  bundle (or a direct sum) on an orthogonal Grassmannian OGr(k, m).  For each
#  one we report:
#    • the full Hodge diamond  h^{p,q}(X);
#    • local rigidity          h¹(X, T_X) = 0
#  Rigidity is read off the LES of  0 → T_X → T_{G/P}|_X → E|_X → 0:
#  when both h¹(T_{G/P}|_X) and h¹(E|_X) vanish (Koszul + Bott give them
#  exactly), the H⁰ surjectivity from quasi-homogeneity yields h¹(T_X) = 0.
#
#  Usage:
#    julia --project=. examples/SpinorialFanos.jl
#    julia --project=. examples/SpinorialFanos.jl --dry-run     # dim/rank only
#    julia --project=. examples/SpinorialFanos.jl --max-dim=14  # cap dim Z
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Printf

# ───────────────────────────────────────────────────────────────────────────────
#  Notation alignment with the paper
#
#  Frassineti–Manivel:  S₊ ↔ ω_n  (highest spin node),  S₋ ↔ ω_{n-1}.
#  Codebase:            spinor_bundle(X, :plus) → ω_{n-1},  :minus → ω_n.
#  Hence F-M's S₊ corresponds to the :minus argument here.  Type B has a
#  unique spinor representation, so the symbol is ignored there.
# ───────────────────────────────────────────────────────────────────────────────

S₊(X) = _is_typeD(X) ? spinor_bundle(X, :minus) : spinor_bundle(X)
S₋(X) = _is_typeD(X) ? spinor_bundle(X, :plus) : spinor_bundle(X)
_is_typeD(X) = dynkin_type(marked_dynkin_type(X)) <: TypeD

# ───────────────────────────────────────────────────────────────────────────────
#  Ambient varieties the codebase's `OGr` does not cover
#
#  OGr(n-1, 2n) = D_n / P_{n-1, n} is the (n-1)-isotropic Grassmannian (Picard
#  rank 2).  The codebase's `OGr(n-1, 2n)` wrongly returns D_n / P_{n-1} (the
#  spinor variety), so we build OGr(3, 8) directly.
#
#  `OGr(n, 2n)` returns D_n/P_n = OGr(n, 2n)₊; the dual component is D_n/P_{n-1}.
# ───────────────────────────────────────────────────────────────────────────────

OGr3_8() = partial_flag_variety(TypeD{4}, (3, 4), "OGr(3, 8)")
OGr_minus(n::Integer) =
  partial_flag_variety(TypeD{Int(n)}, (Int(n) - 1,), "OGr($n, $(2n))_-")

# ───────────────────────────────────────────────────────────────────────────────
#  Table 2
# ───────────────────────────────────────────────────────────────────────────────

struct Variety
  label::String
  build_X::Function          # () -> PartialFlagVariety
  build_E::Function          # X -> CompletelyReducibleBundle
  dim::Int
  picard_rank::Int
  fano_index::Int
  paper_betti::Union{Nothing,Vector{Int}}
end

const TABLE2 = [
  Variety(
    "(OGr(3,8), S₊)", () -> OGr3_8(), X -> S₊(X), 8, 2, 1, [1, 2, 3, 4, 4, 4, 3, 2, 1]
  ),
  Variety(
    "(OGr(3,8), S₊⊕S₋)",
    () -> OGr3_8(),
    X -> S₊(X) ⊕ S₋(X),
    7,
    2,
    3,
    [1, 2, 3, 3, 3, 3, 2, 1],
  ),
  Variety(
    "(OGr(3,9), S)",
    () -> OGr(3, 9),
    X -> spinor_bundle(X),
    10,
    1,
    1,
    [1, 1, 2, 3, 3, 4, 3, 3, 2, 1, 1],
  ),
  Variety(
    "(OGr(2,10), S₊)", () -> OGr(2, 10), X -> S₊(X), 9, 1, 5, [1, 1, 2, 3, 3, 3, 3, 2, 1, 1]
  ),
  Variety(
    "(OGr(3,10), S₊)",
    () -> OGr(3, 10),
    X -> S₊(X),
    13,
    1,
    5,
    [1, 1, 3, 4, 6, 7, 8, 8, 7, 6, 4, 3, 1, 1],
  ),
  Variety("(OGr(3,10), S₊⊕S₊)", () -> OGr(3, 10), X -> S₊(X) ⊕ S₊(X), 11, 1, 4, nothing),
  Variety("(OGr(3,10), S₊⊕S₋)", () -> OGr(3, 10), X -> S₊(X) ⊕ S₋(X), 11, 1, 1, nothing),
  Variety("(OGr(5,10)₊, S₊)", () -> OGr(5, 10), X -> S₊(X), 9, 1, 7, nothing),
  Variety("(OGr(5,10)₊, S₊⊕S₊)", () -> OGr(5, 10), X -> S₊(X) ⊕ S₊(X), 8, 1, 6, nothing),
  Variety("(OGr(3,11), S)", () -> OGr(3, 11), X -> spinor_bundle(X), 14, 1, 6, nothing),
  Variety("(OGr(3,12), S₊)", () -> OGr(3, 12), X -> S₊(X), 17, 1, 7, nothing),
  Variety("(OGr(4,12), S₊)", () -> OGr(4, 12), X -> S₊(X), 20, 1, 6, nothing),
  Variety("(OGr(6,12)₊, S₊)", () -> OGr(6, 12), X -> S₊(X), 14, 1, 9, nothing),
  Variety("(OGr(3,14), S₊)", () -> OGr(3, 14), X -> S₊(X), 19, 1, 6,
    [1, 1, 2, 3, 6, 7, 10, 12, 15, 15, 15, 15, 12, 10, 7, 6, 3, 2, 1, 1]),
  Variety("(OGr(4,14), S₊)", () -> OGr(4, 14), X -> S₊(X), 26, 1, 7, nothing),
  Variety("(OGr(5,14), S₊)", () -> OGr(5, 14), X -> S₊(X), 28, 1, 7, nothing),
  Variety("(OGr(7,14)₊, S₊)", () -> OGr(7, 14), X -> S₊(X), 20, 1, 11, nothing),
  Variety("(OGr(7,14)₋, S₊)", () -> OGr_minus(7), X -> S₊(X), 14, 1, 7, nothing),
]

# ───────────────────────────────────────────────────────────────────────────────
#  Pre-flight: dim G/P and rank E
# ───────────────────────────────────────────────────────────────────────────────

function dry_run()
  rows = Vector{Any}[]
  for v in TABLE2
    X = v.build_X();
    E = v.build_E(X)
    d_amb, r_E = dimension(X), Int(rank_bundle(E))
    push!(
      rows,
      [v.label, d_amb, r_E, d_amb - r_E, v.dim,
        d_amb - r_E == v.dim ? "✓" : "MISMATCH"],
    )
  end
  pretty_table(
    permutedims(reduce(hcat, rows), (2, 1));
    column_labels=["variety", "dim G/P", "rank E", "dim Z", "paper", "status"],
    alignment=[:l, :r, :r, :r, :r, :c],
    display_size=(-1, -1),
  )
end

# ───────────────────────────────────────────────────────────────────────────────
#  Per-variety computation
# ───────────────────────────────────────────────────────────────────────────────

struct Result
  v::Variety
  Z::ZeroLocus
  hodge::Matrix{AffineExpr}
  h1_TX::BigInt              # h¹(Z, T_{G/P}|_Z),  for rigidity
  h1_E::BigInt               # h¹(Z, E|_Z),        for rigidity
  hodge_seconds::Float64
  les_seconds::Float64
end

function compute(v::Variety)
  X = v.build_X()
  E = v.build_E(X)
  Z = zero_locus(E)
  t_h = @elapsed H = hodge_numbers(Z)
  t_l = @elapsed begin
    (H_T, _) = cohomology_on_restriction(Z, tangent_bundle(X))
    (H_E, _) = cohomology_on_restriction(Z, E)
  end
  Result(v, Z, H, H_T[1], H_E[1], t_h, t_l)
end

is_rigid(r::Result) = r.h1_TX == 0 && r.h1_E == 0

# ───────────────────────────────────────────────────────────────────────────────
#  Formatting
# ───────────────────────────────────────────────────────────────────────────────

fmt(e::AffineExpr) = is_determined(e) ? string(Int(e.constant)) : string(e)

diagonal(H::Matrix{AffineExpr}) = [H[p + 1, p + 1] for p in 0:(size(H, 1) - 1)]

function off_diagonal_nonzero(H::Matrix{AffineExpr})
  d = size(H, 1) - 1
  [
    (p, q, H[p + 1, q + 1]) for p in 0:d for q in 0:d
    if p != q && !is_zero_expr(H[p + 1, q + 1])
  ]
end

function free_parameters(H::AbstractMatrix{AffineExpr})
  vars = Set{Int}()
  for e in H, k in keys(e.coeffs)
    push!(vars, k)
  end
  vars
end

# Print the Hodge diamond at any width by using a wide IOContext.
print_diamond_wide(io::IO, H::Matrix{AffineExpr}) =
  print_hodge_diamond(IOContext(io, :displaysize => (10_000, 10_000)), H)

# ───────────────────────────────────────────────────────────────────────────────
#  Reporting
# ───────────────────────────────────────────────────────────────────────────────

function show_one(r::Result)
  v, d = r.v, dimension(r.Z)
  println("─── ", v.label, " ───")
  @printf("  dim Z = %d,  ι_X = %d,  ρ_X = %d   (Hodge %.2fs, LES %.2fs)\n",
    d, v.fano_index, v.picard_rank, r.hodge_seconds, r.les_seconds)

  diag = diagonal(r.hodge)
  println("  h^{p,p} = (", join(fmt.(diag), ", "), ")")
  if v.paper_betti !== nothing
    matches = all(is_determined, diag) &&
              [Int(e.constant) for e in diag] == v.paper_betti
    println(
      matches ? "    ✓ matches paper." :
      "    paper: (" * join(v.paper_betti, ",") * ")",
    )
  end

  off = off_diagonal_nonzero(r.hodge)
  if isempty(off)
    println("  All h^{p,q} = 0 for p ≠ q (cohomology pure).")
  else
    println("  Off-diagonal entries (p ≠ q, not identically zero):")
    for (p, q, e) in off
      println("    h^{$p,$q} = ", fmt(e))
    end
  end

  vars = free_parameters(r.hodge)
  if !isempty(vars)
    n_undet = count(!is_determined, r.hodge)
    @printf("  (%d Hodge entries undetermined, %d free parameter%s)\n",
      n_undet, length(vars), length(vars) == 1 ? "" : "s")
  end

  println("  Hodge diamond:")
  print_diamond_wide(stdout, r.hodge)

  println("  LES rigidity:  h¹(T_{G/P}|_Z) = ", r.h1_TX,
    ",  h¹(E|_Z) = ", r.h1_E)
  println(
    if is_rigid(r)
      "    ⇒ h¹(T_Z) = 0  (locally rigid)"
    else
      "    ⇒ inconclusive from LES dimensions alone"
    end,
  )
  println()
end

function show_summary(results::Vector{Result})
  println("══════════════ Summary ══════════════")
  rows = Vector{Any}[]
  for r in results
    diag = diagonal(r.hodge)
    n_undet = count(!is_determined, r.hodge)
    push!(
      rows,
      [
        r.v.label,
        dimension(r.Z),
        is_rigid(r) ? "yes" : "no",
        all(is_determined, diag) ? "yes" : "no",
        n_undet,
      ],
    )
  end
  pretty_table(
    permutedims(reduce(hcat, rows), (2, 1));
    column_labels=["variety", "dim", "rigid?", "pure diagonal?", "undet. entries"],
    alignment=[:l, :r, :c, :c, :r],
    display_size=(-1, -1),
  )

  n = length(results)
  n_rigid = count(is_rigid, results)
  @printf("\nLocally rigid:  %d / %d\n", n_rigid, n)
end

# ───────────────────────────────────────────────────────────────────────────────
#  Main
# ───────────────────────────────────────────────────────────────────────────────

function main()
  max_dim = typemax(Int)
  for a in ARGS
    m = match(r"^--max-dim=(\d+)$", a)
    m !== nothing && (max_dim = parse(Int, m.captures[1]))
  end

  println("Spinorial Fano manifolds — Frassineti–Manivel arXiv:2605.28712, Table 2\n")
  dry_run()
  println()

  any(==("--dry-run"), ARGS) && return nothing

  results = Result[]
  for v in sort(TABLE2; by=v -> v.dim)
    v.dim > max_dim && continue
    try
      push!(results, compute(v))
      show_one(last(results))
    catch err
      @warn "$(v.label) failed" exception=(err, catch_backtrace())
    end
  end
  show_summary(results)
end

main()
