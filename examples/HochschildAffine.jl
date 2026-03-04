# ═══════════════════════════════════════════════════════════════════════════════
#  HochschildAffine.jl — Polyvector fields on affine cones
#
#  Computes ∧ᵖ T_{G/P} and its cohomology via the Borel–Weil–Bott theorem
#  for a list of generalized Grassmannians G/P.
#
#  Mimics the output of the polyvectors-g-mod-p SageMath script.
#
#  Usage:
#    julia --project=. examples/HochschildAffine.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using Lie
using PrettyTables

# ─── Formatting helpers ─────────────────────────────────────────────────────

"""Format a WeightLatticeElem as ω₁ + 2ω₂ style."""
function format_weight(λ::WeightLatticeElem)
  v = Lie.coefficients(λ)
  parts = String[]
  for (i, c) in enumerate(v)
    c == 0 && continue
    if c == 1
      push!(parts, "ω$i")
    elseif c == -1
      push!(parts, "-ω$i")
    else
      push!(parts, "$(c)ω$i")
    end
  end
  isempty(parts) && return "0"
  s = join(parts, " + ")
  s = replace(s, "+ -" => "- ")
  s
end

"""Superscript digits for display."""
function superscript(n::Int)
  digits = Dict(
    '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
  )
  String([get(digits, c, c) for c in string(n)])
end

# ─── Main computation ────────────────────────────────────────────────────────

function polyvectors_table(X::PartialFlagVariety)
  MDT = marked_type(X)
  DT = dynkin_type(X)
  d = dimension(X)

  T = tangent_bundle(X)
  n_bundles = n_components(T)

  println("For $X (tangent bundle decomposes into $n_bundles summand(s))")
  println()

  total_euler = BigInt(0)

  for p in 0:d
    Ep = exterior_power(T, p)
    println("∧$(superscript(p)) T  (rank = $(rank_bundle(Ep)), expected = $(binomial(d, p)))")

    # Collect rows for the table
    weights = String[]
    ranks = BigInt[]
    regulars = String[]
    dominants = String[]
    degrees = String[]
    dims = String[]

    for comp in components(Ep)
      λ = to_ambient_weight(MDT, comp)
      r = fiber_dimension(comp)

      push!(weights, format_weight(λ))
      push!(ranks, r)

      bwb = borel_weil_bott(λ)
      if bwb === nothing
        push!(regulars, "false")
        push!(dominants, "")
        push!(degrees, "")
        push!(dims, "")
      else
        deg, μ = bwb
        push!(regulars, "true")
        push!(dominants, format_weight(μ))
        push!(degrees, string(deg))
        push!(dims, string(degree(μ)))
      end
    end

    data = hcat(weights, ranks, regulars, dominants, degrees, dims)
    pretty_table(data;
      column_labels=["weight", "rank", "regular", "dominant", "degree", "dimension"],
      alignment=[:l, :r, :c, :l, :r, :r],
    )

    n_summands = length(components(Ep))
    println("  $n_summands summand(s)")

    # Compute Euler characteristic for this exterior power
    H = dimensions(Ep)
    χ = euler_characteristic(H)
    total_euler += (-1)^p * χ

    println()
  end

  # Summary
  χ_X = euler_characteristic(X)
  println("─── Summary for $X ───")
  println("  Dimension of G/P: $d")
  println("  Rank of K₀:       $χ_X")
  expected = (-1)^d * χ_X
  println("  Litmus test:      $total_euler == $expected, $(total_euler == expected)")
  println()
end

# ─── Run on standard examples ───────────────────────────────────────────────

function main()
  varieties = [
    # Type A
    partial_flag_variety(TypeA{2}, 1),
    partial_flag_variety(TypeA{3}, 2),
    # Type B
    partial_flag_variety(TypeB{2}, 1),
    # Type G₂
    partial_flag_variety(TypeG2, 1),
    partial_flag_variety(TypeG2, 2),
    # Type E₆
    partial_flag_variety(TypeE{6}, 1),
  ]

  for X in varieties
    println("=" ^ 72)
    polyvectors_table(X)
  end
end

main()
