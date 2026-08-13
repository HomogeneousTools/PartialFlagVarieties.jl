# ═══════════════════════════════════════════════════════════════════════════════
#  Kuechle.jl — Fano fourfolds classified by Küchle and Fatighenti–Mongardi
#
#  Reproduces the numerical invariants from:
#    Küchle, "On Fano 4-folds of index 1 and homogeneous vector bundles
#    over Grassmannians", Math. Z. 218 (1995), 563–575.
#  and the additional K3-type families from:
#    Fatighenti–Mongardi, "Fano varieties of K3 type and IHS manifolds",
#    Int. Math. Res. Not. (2021), arXiv:1904.05679.
#
#  Küchle's 21 Fano-index-1 families are labeled b1–d3 in his notation.
#  h^0(-K_Z) is computed via sheaf cohomology: by adjunction
#    ω_Z ≅ (ω_X ⊗ det(E))|_Z,
#  so H^0(Z, ω_Z^{-1}) = H^0(Z, (anticanonical(X) ⊗ det(E)^{-1})|_Z),
#  computed through the Koszul long exact sequence.
#  Hodge numbers come from the same Koszul long exact sequence.
#
#  Usage: julia --project=. examples/Kuechle.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Printf

# =============================================================================
#  GL(n) weight -> omega-basis and bundle construction
# =============================================================================

"Convert GL(n) weight for Gr(k,n) to fundamental weight basis (Q first)."
function gl_weight_to_omega(k::Int, n::Int, w::Vector{Int})
  q_part = w[1:(n - k)]
  s_part = w[(n - k + 1):n]
  eps = vcat(s_part, q_part)
  [eps[i] - eps[i + 1] for i in 1:(n - 1)]
end

"Build CompletelyReducibleBundle on Gr(k,n) from a list of GL(n) weight vectors."
function bundle_from_gl_weights(
  X::PartialFlagVariety, k::Int, n::Int, weights::Vector{Vector{Int}}
)
  summands = Vector{Int}[]
  for w in weights
    push!(summands, gl_weight_to_omega(k, n, w))
  end
  CompletelyReducibleBundle(X, summands)
end

# =============================================================================
#  Fano index
# =============================================================================

# Fano index of Z(E) in G/P via adjunction, or -1 if Picard rank > 1.
function _fano_index_zl(Z)
  try
    fano_index(Z)
  catch ArgumentError
    -1
  end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Fano fourfold families in Küchle's notation (Theorem 3.1)
#
#  Format: (label, k, n, [GL(n)-weight-vectors], description)
#
#  GL(n) weights: first (n-k) entries = Q-side, last k entries = S-side.
#  Omitted: c4, d1, d2 (not of interest).
# ═══════════════════════════════════════════════════════════════════════════════

# Küchle's 21 Fano-index-1 families, labeled b1–d3
const KUECHLE_FAMILIES = [
  # ── Gr(2,5): b1, b2 ────────────────────────────────────────────────
  ("b1", 2, 5, [[0, 0, 0, 3, 3], [0, 0, 0, 1, 1]], "O(3) + O(1) on Gr(2,5)"),
  ("b2", 2, 5, [[0, 0, 0, 2, 2], [0, 0, 0, 2, 2]], "O(2)² on Gr(2,5)"),

  # ── Gr(2,6): b3–b6 ─────────────────────────────────────────────────
  ("b3", 2, 6, [[1, 1, 1, 0, 2, 2]], "∧³Q* ⊗ O(2) on Gr(2,6)"),
  ("b4", 2, 6, [[0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 2, 2]], "Sym²S* + O(2) on Gr(2,6)"),
  ("b5", 2, 6, [[0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 2, 1]],
    "O(1)² + S*(1) on Gr(2,6)"),
  ("b6", 2, 6,
    [[0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 2, 2]],
    "O(1)³ + O(2) on Gr(2,6)"),

  # ── Gr(2,7): b7–b10 ────────────────────────────────────────────────
  ("b7", 2, 7,
    [[0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 1],
      [0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 1]],
    "O(1)⁶ on Gr(2,7)"),
  ("b8", 2, 7,
    [
      [0, 0, 0, 0, 0, 1, 1],
      [0, 0, 0, 0, 0, 1, 1],
      [0, 0, 0, 0, 0, 1, 1],
      [0, 0, 0, 0, 0, 2, 0],
    ],
    "O(1)³ + Sym²S* on Gr(2,7)"),
  ("b9", 2, 7, [[0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 2, 0]], "(Sym²S*)² on Gr(2,7)"),
  ("b10", 2, 7, [[0, 0, 0, 0, 0, 2, 2], [1, 0, 0, 0, 0, 1, 1]], "O(2) + Q*(1) on Gr(2,7)"),

  # ── Gr(2,8): b11 ───────────────────────────────────────────────────
  ("b11", 2, 8,
    [[0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 1, 1], [1, 0, 0, 0, 0, 0, 1, 1]],
    "O(1)² + Q*(1) on Gr(2,8)"),

  # ── Gr(3,6): c1, c2 ────────────────────────────────────────────────
  ("c1", 3, 6,
    [[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1],
      [0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
    "O(1)⁵ on Gr(3,6)"),
  ("c2", 3, 6,
    [[0, 0, 0, 1, 1, 0], [0, 0, 0, 1, 1, 1], [0, 0, 0, 2, 2, 2]],
    "∧²S* + O(1) + O(2) on Gr(3,6)"),

  # ── Gr(3,7): c3–c6 ─────────────────────────────────────────────────
  ("c3", 3, 7, [[1, 0, 0, 0, 1, 1, 1], [1, 0, 0, 0, 1, 1, 1]], "Q*(1)² on Gr(3,7)"),
  ("c5", 3, 7,
    [[0, 0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 1, 1, 1], [1, 0, 0, 0, 1, 1, 1]],
    "∧²S* + O(1) + Q*(1) on Gr(3,7)"),
  ("c6", 3, 7,
    [
      [0, 0, 0, 0, 1, 1, 0],
      [0, 0, 0, 0, 1, 1, 0],
      [0, 0, 0, 0, 1, 1, 1],
      [0, 0, 0, 0, 1, 1, 1],
    ],
    "(∧²S*)² + O(1)² on Gr(3,7)"),

  # ── Gr(3,8): c7 ────────────────────────────────────────────────────
  ("c7", 3, 8,
    [[0, 0, 0, 0, 0, 1, 1, 1], [1, 1, 0, 0, 0, 1, 1, 1]],
    "O(1) + (∧²Q* ⊗ O(1)) on Gr(3,8)"),

  # ── Gr(5,10): d3 ───────────────────────────────────────────────────
  #("d3", 5, 10,
  #  [[0, 0, 0, 0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]],
  #  "(∧²S*)² + O(1) on Gr(5,10)"),
]

# Fatighenti–Mongardi K3-type families (Fano index > 1, not in Küchle)
const FM_FAMILIES = [
  ("FM-i", 4, 8,
    [[0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 0, 0]],
    "(∧²S*)² on Gr(4,8)"),
  ("FM-ii", 4, 8,
    [[0, 0, 0, 0, 1, 1, 0, 0], [1, 1, 0, 0, 1, 1, 1, 1]],
    "∧²S* + (∧²Q* ⊗ O(1)) on Gr(4,8)"),
  ("FM-iii", 5, 10,
    [
      [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
      [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
      [1, 1, 1, 0, 0, 1, 1, 1, 1, 1],
    ],
    "∧²S* + O(1) + (∧³Q* ⊗ O(1)) on Gr(5,10)"),
]

# ═══════════════════════════════════════════════════════════════════════════════
#  Result type
# ═══════════════════════════════════════════════════════════════════════════════

struct FanoResult
  label::String
  description::String
  k::Int
  n::Int
  h0_antiK::Union{BigInt,Nothing}
  chi_top::Any
  b2::Any
  b3::Any
  h11::Any
  h22::Any
  h13::Any
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Per-family computation
# ═══════════════════════════════════════════════════════════════════════════════

function compute_family(label, k, n, weights, desc)
  X = Gr(k, n)
  d = dimension(X)
  E = bundle_from_gl_weights(X, k, n, weights)
  r = Int(rank(E))

  if d - r != 4
    @warn "  $label: dim(Z) = $(d - r) ≠ 4; skipping"
    return nothing
  end

  Z = zero_locus(E)

  # h^0(-K_Z) via sheaf cohomology: by adjunction ω_Z ≅ (ω_X ⊗ det(E))|_Z,
  # so H^0(Z, ω_Z^{-1}) = H^0(Z, F_antiK|_Z) where
  #   F_antiK = anticanonical_bundle(X) ⊗ dual(det(E)).
  h0aK = nothing
  try
    F_antiK = anticanonical_bundle(X) ⊗ dual(det(E))
    (H_antiK, _) = cohomology_on_restriction(Z, F_antiK)
    h0aK = H_antiK[0]
  catch e
    @warn "  $label h⁰(-K) failed: $(sprint(showerror, e))"
  end

  # Hodge numbers, Betti numbers, topological Euler characteristic
  chi = nothing
  β2 = nothing
  β3 = nothing
  h11 = nothing
  h22 = nothing
  h13 = nothing
  try
    hm = hodge_numbers(Z)
    # χ_top = Σ_{p,q} (-1)^{p+q} h^{p,q}
    chi = sum((-1)^(p + q) * hm[p + 1, q + 1] for p in 0:4, q in 0:4)
    β2 = hm[2, 2]                    # h^{1,1}
    β3 = hm[3, 2] + hm[2, 3]        # h^{2,1} + h^{1,2} = 2 h^{2,1}
    h11 = hm[2, 2]                   # h^{1,1}
    h22 = hm[3, 3]                   # h^{2,2}
    h13 = hm[2, 4]                   # h^{1,3}
  catch e
    @warn "  $label Hodge failed: $(sprint(showerror, e))"
  end

  FanoResult(label, desc, k, n, h0aK, chi, β2, β3, h11, h22, h13)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Display helpers
# ═══════════════════════════════════════════════════════════════════════════════

fmt(::Nothing) = "—"
fmt(x) = string(x)

function show_kuechle_table(results)
  rows = map(results) do r
    [r.label, "Gr($(r.k),$(r.n))", r.description,
      fmt(r.h0_antiK), fmt(r.chi_top), fmt(r.b2), fmt(r.b3),
      fmt(r.h11), fmt(r.h22), fmt(r.h13)]
  end
  data = permutedims(hcat(rows...), (2, 1))
  pretty_table(
    data;
    column_labels=[
      "#", "G/P", "Bundle E", "h⁰(-K)", "χ_top", "b₂", "b₃", "h¹¹", "h²²", "h¹³"
    ],
    alignment=[:l, :c, :l, :r, :r, :r, :r, :r, :r, :r],
    display_size=(-1, -1),
  )
end

function show_fm_table(results)
  println()
  rows = map(results) do r
    [r.label, "Gr($(r.k),$(r.n))", r.description, string(r.fano_index),
      fmt(r.h0_antiK), fmt(r.chi_top), fmt(r.b2), fmt(r.b3),
      fmt(r.h11), fmt(r.h22), fmt(r.h13)]
  end
  data = permutedims(hcat(rows...), (2, 1))
  pretty_table(
    data;
    column_labels=[
      "#", "G/P", "Bundle E", "h⁰(-K)", "χ_top", "b₂", "b₃", "h¹¹", "h²²", "h¹³"
    ],
    alignment=[:l, :c, :l, :r, :r, :r, :r, :r, :r, :r, :r],
    display_size=(-1, -1),
  )
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

function main()
  println("\nComputing Küchle fourfold invariants...\n")

  kuechle_results = FanoResult[]
  t_total = @elapsed for (label, k, n, weights, desc) in KUECHLE_FAMILIES
    print("  $label ... ")
    flush(stdout)
    res = compute_family(label, k, n, weights, desc)
    if res !== nothing
      push!(kuechle_results, res)
      println("done")
    end
  end
  println(@sprintf("\n  Total: %.2fs", t_total))

  #println("\nComputing Fatighenti–Mongardi families...")
  #fm_results = FanoResult[]
  #for (label, k, n, weights, desc) in FM_FAMILIES
  #  print("  $label ... ")
  #  flush(stdout)
  #  res = compute_family(label, k, n, weights, desc)
  #  if res !== nothing
  #    push!(fm_results, res)
  #    println("done")
  #  end
  #end

  println()
  show_kuechle_table(kuechle_results)
  #show_fm_table(fm_results)
end

main()
