# ═══════════════════════════════════════════════════════════════════════════════
#  BottVanishing.jl — Bott vanishing on (co)adjoint varieties
#
#  For each (co)adjoint variety G/P with Picard rank 1, computes the
#  smallest positive q such that H^p(G/P, Ω^i(q)) ≠ 0 for some p > 0.
#
#  A smooth projective variety X satisfies Bott vanishing if
#    H^p(X, Ω^i_X ⊗ L) = 0  for all p > 0, i ≥ 0, and ample L.
#  This script finds the first q where this fails.
#
#  Related to arXiv:2506.09811.
#
#  Usage:
#    julia --project=. examples/BottVanishing.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using Semisimple
using PrettyTables

"""
    find_bott_vanishing_q(X; max_q=20) -> Int

For a partial flag variety X with Picard rank 1, find the smallest q ≥ 1
such that H^p(X, Ω^i(q)) ≠ 0 for some p > 0 and some 0 ≤ i ≤ dim(X).

Returns the q, or -1 if no failure found up to max_q.
"""
function find_bott_vanishing_q(X; max_q=20)
  d = dimension(X)
  Ω = cotangent_bundle(X)

  for q in 1:max_q
    for i in 0:d
      Ei = exterior_power(Ω, i)
      Et = twist(Ei, 1, q)
      H = dimensions(Et)

      for p in 1:d
        if H[p] != 0
          return q
        end
      end
    end
  end
  return -1
end

"""
    find_bott_vanishing_witnesses(X, q) -> Vector

For the given q, find all (i, p) such that H^p(X, Ω^i(q)) ≠ 0 and p > 0.
"""
function find_bott_vanishing_witnesses(X, q)
  d = dimension(X)
  Ω = cotangent_bundle(X)
  witnesses = Tuple{Int,Int,BigInt}[]

  for i in 0:d
    Ei = exterior_power(Ω, i)
    Et = twist(Ei, 1, q)
    H = dimensions(Et)
    for p in 1:d
      if H[p] != 0
        push!(witnesses, (i, p, H[p]))
      end
    end
  end
  return witnesses
end

function main()
  # Adjoint varieties with Picard rank 1
  # Type A adjoint varieties have Picard rank 2, so are excluded
  adjoint_cases = [
    ("G₂", adjoint_variety(TypeG2)),
    ("B₃", adjoint_variety(TypeB{3})),
    ("D₄", adjoint_variety(TypeD{4})),
    ("B₄", adjoint_variety(TypeB{4})),
    ("D₅", adjoint_variety(TypeD{5})),
    ("F₄", adjoint_variety(TypeF4)),
    ("E₆", adjoint_variety(TypeE{6})),
    # E₇ is expensive (dim 33), uncomment to run:
    # ("E₇", adjoint_variety(TypeE{7})),
    # E₈ is infeasible: ∧^p of a rank-56 bundle
  ]

  println("Bott vanishing failure for adjoint varieties")
  println("=" ^ 60)
  println()

  labels = String[]
  varieties_str = String[]
  dims = Int[]
  q_vals = Int[]

  for (name, X) in adjoint_cases
    d = dimension(X)
    print("Computing $name (dim=$d)... ")
    flush(stdout)

    q = find_bott_vanishing_q(X; max_q=15)
    push!(labels, name)
    push!(varieties_str, string(X))
    push!(dims, d)
    push!(q_vals, q)

    println("q = $q")

    # Show witnesses
    ws = find_bott_vanishing_witnesses(X, q)
    for (i, p, dim_val) in ws
      println("    H$p(Ω^$i($q)) = $dim_val")
    end
  end

  println()
  data = hcat(labels, varieties_str, dims, q_vals)
  pretty_table(data;
    column_labels=["Type", "Variety", "dim", "q (failure)"],
    alignment=[:l, :l, :r, :r],
  )

  # Coadjoint varieties (where different from adjoint)
  println()
  println("Bott vanishing failure for coadjoint varieties")
  println("=" ^ 60)
  println()

  coadjoint_cases = [
    ("C₃", coadjoint_variety(TypeC{3})),
    ("C₄", coadjoint_variety(TypeC{4})),
    ("F₄", coadjoint_variety(TypeF4)),
  ]

  for (name, X) in coadjoint_cases
    d = dimension(X)
    print("Computing $name (dim=$d)... ")
    flush(stdout)

    q = find_bott_vanishing_q(X; max_q=5)
    println("q = $q")

    ws = find_bott_vanishing_witnesses(X, q)
    for (i, p, dim_val) in ws
      println("    H$p(Ω^$i($q)) = $dim_val")
    end
  end
end

main()
