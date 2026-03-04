# ═══════════════════════════════════════════════════════════════════════════════
#  Kuechle.jl — Fano fourfolds as zero loci in Grassmannians
#
#  Classifies Fano fourfolds arising as zero loci of sections of
#  equivariant vector bundles on Grassmannians Gr(k,n).
#
#  Each family is specified by a list of GL(n) weight vectors that encode
#  the irreducible summands of the defining bundle E.  A GL(n) weight
#  vector (d₁,…,dₙ) for Gr(k,n) has:
#    • Q-side: first n−k entries (quotient bundle part)
#    • S-side: last k entries (tautological subbundle part)
#  and converts to an ω-basis weight for A_{n−1} via the change of basis
#    ε = (S-side, Q-side),   ω_i = ε_i − ε_{i+1}.
#
#  Reference: Küchle, "On Fano 4-folds of index 1 and homogeneous vector
#             bundles over Grassmannians" (Math. Z. 218, 1995)
#  Extended:  Fatighenti–Mongardi, "Fano varieties of K3 type and IHS
#             manifolds" (arXiv:1904.05679)
#
#  Usage: julia --project=. examples/Kuechle.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Lie

# ═══════════════════════════════════════════════════════════════════════════════
#  GL(n) weight → ω-basis weight conversion
# ═══════════════════════════════════════════════════════════════════════════════

"""
Convert a GL(n) weight vector `w` for Gr(k,n) to the fundamental weight
basis of A_{n-1}.

The GL(n) weight is split as (Q-side, S-side) = (w[1:n-k], w[n-k+1:n]).
The ε-vector is (S-side, Q-side) in the standard GL(n) ordering:
  ε₁,…,εₖ correspond to S, εₖ₊₁,…,εₙ correspond to Q.
Then ωᵢ = εᵢ − εᵢ₊₁ for i = 1,…,n−1.
"""
function gl_weight_to_omega(k::Int, n::Int, w::Vector{Int})
  @assert length(w) == n
  q_part = w[1:n-k]
  s_part = w[n-k+1:n]
  eps = vcat(s_part, q_part)  # S-side first, then Q-side
  [eps[i] - eps[i+1] for i in 1:n-1]
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Bundle construction from GL(n) weights
# ═══════════════════════════════════════════════════════════════════════════════

"""
Construct a CompletelyReducibleBundle on Gr(k,n) from a list of GL(n)
weight vectors.  Each weight vector specifies an irreducible summand.
"""
function bundle_from_gl_weights(X::PartialFlagVariety{MDT}, k::Int, n::Int, weights::Vector{Vector{Int}}) where {MDT}
  DT = PartialFlagVarieties._ambient_type(MDT)

  summands = IrrepLevi{MDT}[]
  for w in weights
    omega = gl_weight_to_omega(k, n, w)
    λ = WeightLatticeElem(DT, omega)
    push!(summands, IrrepLevi(MDT, λ))
  end

  CompletelyReducibleBundle{MDT}(X, summands)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Fano fourfold families — data from Küchle (1995) and extensions
#
#  Format: (label, k, n, [weight_vectors], description)
#
#  Labels follow the convention:
#    a = hypersurfaces and complete intersections in projective spaces
#    b = bundles involving S* (tautological dual)
#    c = bundles involving Q (quotient) or mixed
#    d = larger Grassmannians
# ═══════════════════════════════════════════════════════════════════════════════

const FAMILIES = [
  # ── Hypersurfaces in P^5 ────────────────────────────────────────────
  ("a1",  1, 6, [[0,0,0,0,0,2]], "O(2) on ℙ⁵"),
  ("a2",  1, 6, [[0,0,0,0,0,3]], "O(3) on ℙ⁵"),
  ("a3",  1, 6, [[0,0,0,0,0,4]], "O(4) on ℙ⁵"),
  ("a4",  1, 6, [[0,0,0,0,0,5]], "O(5) on ℙ⁵"),

  # ── Complete intersections in P^6 ───────────────────────────────────
  ("a5",  1, 7, [[0,0,0,0,0,0,2], [0,0,0,0,0,0,2]], "(2,2) CI in ℙ⁶"),
  ("a6",  1, 7, [[0,0,0,0,0,0,2], [0,0,0,0,0,0,3]], "(2,3) CI in ℙ⁶"),
  ("a7",  1, 7, [[0,0,0,0,0,0,2], [0,0,0,0,0,0,4]], "(2,4) CI in ℙ⁶"),
  ("a8",  1, 7, [[0,0,0,0,0,0,3], [0,0,0,0,0,0,3]], "(3,3) CI in ℙ⁶"),

  # ── Gr(2,5) families ────────────────────────────────────────────────
  ("b1",  2, 5, [[0,0,0,2,1]], "S*(1) on Gr(2,5)"),
  ("b2",  2, 5, [[0,0,0,1,1], [0,0,0,1,1]], "O(1)² on Gr(2,5)"),
  ("b3",  2, 5, [[0,0,0,1,1], [0,0,0,2,2]], "O(1)⊕O(2) on Gr(2,5)"),
  ("b4",  2, 5, [[0,0,0,1,1], [0,0,0,3,3]], "O(1)⊕O(3) on Gr(2,5)"),
  ("b5",  2, 5, [[0,0,0,2,2], [0,0,0,2,2]], "O(2)² on Gr(2,5)"),

  # ── Complete intersections in P^7 ───────────────────────────────────
  ("a9",  1, 8, [[0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,2]],
   "(2,2,2) CI in ℙ⁷"),
  ("a10", 1, 8, [[0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,3]],
   "(2,2,3) CI in ℙ⁷"),

  # ── P^8 complete intersection ───────────────────────────────────────
  ("a11", 1, 9,
   [[0,0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,0,2],
    [0,0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,0,2]],
   "(2,2,2,2) CI in ℙ⁸"),

  # ── Gr(2,6) families ────────────────────────────────────────────────
  ("c1",  2, 6, [[1,0,0,0,1,1]], "Q* on Gr(2,6)"),
  ("c2",  2, 6, [[1,1,1,0,2,2]], "∧³Q*⊗O(1) on Gr(2,6)"),
  ("c3",  2, 6, [[0,0,0,0,1,1], [0,0,0,0,2,0]],
   "O(1)⊕Sym²S* on Gr(2,6)"),
  ("c4",  2, 6, [[0,0,0,0,2,0], [0,0,0,0,2,2]],
   "Sym²S*⊕O(2) on Gr(2,6)"),
  ("c5",  2, 6, [[0,0,0,0,1,1], [0,0,0,0,1,1], [0,0,0,0,2,1]],
   "O(1)²⊕S*(1) on Gr(2,6)"),
  ("c6",  2, 6, [[0,0,0,0,1,1], [0,0,0,0,1,1], [0,0,0,0,1,1], [0,0,0,0,1,1]],
   "O(1)⁴ on Gr(2,6)"),
  ("c7",  2, 6, [[0,0,0,0,1,1], [0,0,0,0,1,1], [0,0,0,0,1,1], [0,0,0,0,2,2]],
   "O(1)³⊕O(2) on Gr(2,6)"),

  # ── Gr(3,6) families ────────────────────────────────────────────────
  ("d1",  3, 6, [[0,0,0,1,1,0], [0,0,0,1,1,1], [0,0,0,1,1,1]],
   "∧²S*⊕O(1)² on Gr(3,6)"),
  ("d2",  3, 6, [[0,0,0,1,1,0], [0,0,0,1,1,1], [0,0,0,2,2,2]],
   "∧²S*⊕O(1)⊕O(2) on Gr(3,6)"),
  ("d3",  3, 6,
   [[0,0,0,1,1,1], [0,0,0,1,1,1], [0,0,0,1,1,1], [0,0,0,1,1,1], [0,0,0,1,1,1]],
   "O(1)⁵ on Gr(3,6)"),

  # ── Gr(2,7) families ────────────────────────────────────────────────
  ("e1",  2, 7, [[0,0,0,0,0,1,1], [1,0,0,0,0,1,1]],
   "O(1)⊕Q* on Gr(2,7)"),
  ("e2",  2, 7, [[1,0,0,0,0,1,1], [0,0,0,0,0,2,2]],
   "Q*⊕O(2) on Gr(2,7)"),
  ("e3",  2, 7, [[0,0,0,0,0,2,0], [0,0,0,0,0,2,0]],
   "Sym²S*⊕Sym²S* on Gr(2,7)"),
  ("e4",  2, 7,
   [[0,0,0,0,0,1,1], [0,0,0,0,0,1,1], [0,0,0,0,0,1,1], [0,0,0,0,0,2,0]],
   "O(1)³⊕Sym²S* on Gr(2,7)"),
  ("e5",  2, 7,
   [[0,0,0,0,0,1,1], [0,0,0,0,0,1,1], [0,0,0,0,0,1,1],
    [0,0,0,0,0,1,1], [0,0,0,0,0,1,1], [0,0,0,0,0,1,1]],
   "O(1)⁶ on Gr(2,7)"),

  # ── Gr(2,8) ─────────────────────────────────────────────────────────
  ("f1",  2, 8,
   [[0,0,0,0,0,0,1,1], [0,0,0,0,0,0,1,1], [1,0,0,0,0,0,1,1]],
   "O(1)²⊕Q* on Gr(2,8)"),

  # ── Gr(3,7) families ────────────────────────────────────────────────
  ("g1",  3, 7, [[1,0,0,0,1,1,1], [1,0,0,0,1,1,1]],
   "Q*⊕Q* on Gr(3,7)"),
  ("g2",  3, 7, [[0,0,0,0,1,1,0], [0,0,0,0,1,1,1], [1,0,0,0,1,1,1]],
   "∧²S*⊕O(1)⊕Q* on Gr(3,7)"),
  ("g3",  3, 7, [[0,0,0,0,1,1,1], [0,0,0,0,1,1,1], [0,0,0,0,2,0,0]],
   "O(1)²⊕Sym²S* on Gr(3,7)"),
  ("g4",  3, 7,
   [[0,0,0,0,1,1,0], [0,0,0,0,1,1,0], [0,0,0,0,1,1,1], [0,0,0,0,1,1,1]],
   "∧²S*⊕∧²S*⊕O(1)² on Gr(3,7)"),

  # ── Gr(3,8) ─────────────────────────────────────────────────────────
  ("h1",  3, 8, [[0,0,0,0,0,1,1,1], [1,1,0,0,0,1,1,1]],
   "O(1)⊕Σ²¹Q* on Gr(3,8)"),

  # ── Gr(4,8) families ────────────────────────────────────────────────
  ("i1",  4, 8, [[0,0,0,0,1,1,0,0], [0,0,0,0,1,1,0,0]],
   "∧²S*⊕∧²S* on Gr(4,8)"),
  ("i2",  4, 8, [[0,0,0,0,1,1,0,0], [1,1,0,0,1,1,1,1]],
   "∧²S*⊕∧²Q*⊗O(1) on Gr(4,8)"),
  ("i3",  4, 8,
   [[0,0,0,0,1,1,1,1], [0,0,0,0,1,1,1,1], [0,0,0,0,2,0,0,0]],
   "O(1)²⊕Sym²S* on Gr(4,8)"),

  # ── Gr(4,9) ─────────────────────────────────────────────────────────
  ("j1",  4, 9, [[0,0,0,0,0,1,1,0,0], [0,0,0,0,0,2,0,0,0]],
   "∧²S*⊕Sym²S* on Gr(4,9)"),

  # ── Gr(5,10) families ───────────────────────────────────────────────
  ("k1",  5, 10,
   [[0,0,0,0,0,1,1,0,0,0], [0,0,0,0,0,1,1,0,0,0], [0,0,0,0,0,1,1,1,1,1]],
   "∧²S*⊕∧²S*⊕O(1) on Gr(5,10)"),
  ("k2",  5, 10,
   [[0,0,0,0,0,1,1,0,0,0], [0,0,0,0,0,1,1,1,1,1], [1,1,1,0,0,1,1,1,1,1]],
   "∧²S*⊕O(1)⊕∧³Q*⊗O(1) on Gr(5,10)"),
]

# ═══════════════════════════════════════════════════════════════════════════════
#  Fano index computation
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compute the Fano index of the zero locus Z = Z(E) ⊂ G/P, i.e., the
largest integer r such that -K_Z = r·H for H the positive generator
of Pic(Z).  For Picard rank 1:
  index(Z) = (anticK_central − det_central) / O(1)_central
where all quantities are in central character coordinates.
"""
function fano_index_zero_locus(X, E)
  MDT = marked_type(X)
  Marked = marked_nodes(MDT)
  length(Marked) == 1 || return -1

  m = Marked[1]
  anticK_c = PartialFlagVarieties._anticanonical_central(MDT)
  det_c = PartialFlagVarieties._determinant_central(E)
  M = PartialFlagVarieties.decomposition_matrix(MDT)
  o1_c = M[m, m]  # central character of O(1)

  idx = (anticK_c[1] - det_c[1]) / o1_c
  Int(idx)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Result type
# ═══════════════════════════════════════════════════════════════════════════════

struct FanoResult
  label::String
  description::String
  k::Int
  n::Int
  dim_ambient::Int
  rank_E::BigInt
  fano_index::Int
  chi_O::BigInt
  h11::BigInt
  h22::BigInt
  determined::Bool
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main computation
# ═══════════════════════════════════════════════════════════════════════════════

function main()
  println("=" ^ 72)
  println("  Fano fourfolds as zero loci in Grassmannians")
  println("  Küchle (1995), Fatighenti–Mongardi (2019)")
  println("=" ^ 72)
  println()

  # Maximum codimension for which we attempt Hodge number computation.
  # Higher codimension means more exterior powers and more expensive
  # BWB computations.
  MAX_CODIM_HODGE = 8

  results = FanoResult[]

  for (label, k, n, weights, desc) in FAMILIES
    print("  $label: $desc ... ")
    flush(stdout)

    X = Gr(k, n)
    d = dimension(X)

    # Build the equivariant bundle E
    E = bundle_from_gl_weights(X, k, n, weights)
    r = rank_bundle(E)
    codim = Int(r)

    # Sanity check: zero locus should be a 4-fold
    if d - codim != 4
      println("SKIP (dim Z = $(d - codim) ≠ 4)")
      continue
    end

    Z = zero_locus(E)

    # Fano index of the zero locus
    idx = fano_index_zero_locus(X, E)

    # Euler characteristic χ(O_Z) — always exact via alternating Koszul sum
    chi_O = euler_characteristic(Z)

    # Hodge numbers — attempt only for small codimension
    h11 = BigInt(-1)
    h22 = BigInt(-1)
    det = false

    if codim <= MAX_CODIM_HODGE
      try
        h = hodge_numbers(Z)
        h11 = h[2, 2]  # h^{1,1}
        h22 = h[3, 3]  # h^{2,2}
        det = true
      catch e
        println("(Hodge failed: $(sprint(showerror, e))) ")
      end
    end

    push!(results, FanoResult(label, desc, k, n, d, r, idx,
      chi_O, h11, h22, det))
    println("done (χ(O) = $chi_O)")
  end

  # ── Display results ───────────────────────────────────────────────────
  println()
  println("=" ^ 72)
  println("  Results: $(length(results)) Fano fourfold families")
  println("=" ^ 72)
  println()

  labels = String[]
  descs = String[]
  grassmannians = String[]
  ranks = String[]
  indices = String[]
  chis = String[]
  h11s = String[]
  h22s = String[]

  for r in results
    push!(labels, r.label)
    push!(descs, r.description)
    push!(grassmannians, "Gr($(r.k),$(r.n))")
    push!(ranks, string(r.rank_E))
    push!(indices, string(r.fano_index))
    push!(chis, string(r.chi_O))
    push!(h11s, r.h11 >= 0 ? string(r.h11) : "—")
    push!(h22s, r.h22 >= 0 ? string(r.h22) : "—")
  end

  data = hcat(labels, descs, grassmannians, ranks, indices, chis, h11s, h22s)
  pretty_table(data;
    column_labels=["#", "Bundle", "G/P", "rk(E)", "i(Z)", "χ(O)", "h¹¹", "h²²"],
    alignment=[:l, :l, :c, :r, :r, :r, :r, :r],
    fit_table_in_display_vertically=false,
  )
end

main()
