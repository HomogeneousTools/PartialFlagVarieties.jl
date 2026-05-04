# ═══════════════════════════════════════════════════════════════════════════════
#  HochschildKuechle.jl — Hochschild cohomology of Küchle fourfolds
#
#  Reuses the Küchle families from examples/Kuechle.jl and prints the HKR
#  polyvector parallelogram of each zero locus together with its label and
#  geometric description.
#
#  Usage:
#    julia --project=. examples/HochschildKuechle.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Printf

include(joinpath(@__DIR__, "Kuechle.jl"))

struct HochschildKuechleResult
  label::String
  description::String
  ambient::String
  bundle::String
  hh_dims::Vector{Any}
  polyvector::PolyvectorParallelogram
end

_fmt(x) = string(x)

function _hh_dims(P::PolyvectorParallelogram)
  [hochschild_dimension(P, n) for n in 0:(2 * P.dim)]
end

function compute_hochschild_family(label, k, n, weights, desc)
  X = Gr(k, n)
  E = bundle_from_gl_weights(X, k, n, weights)
  Z = zero_locus(E)
  P = hochschild_cohomology(Z)
  HochschildKuechleResult(label, desc, string(X), desc, _hh_dims(P), P)
end

function show_summary(results)
  rows = map(results) do r
    [r.label, r.ambient, r.description, join(_fmt.(r.hh_dims), ", ")]
  end
  data = permutedims(hcat(rows...), (2, 1))
  pretty_table(
    data;
    column_labels=["#", "G/P", "Bundle E", "HHⁿ dims for n = 0,…,8"],
    alignment=[:l, :c, :l, :l],
    display_size=(-1, -1),
  )
end

function show_result(io::IO, r::HochschildKuechleResult)
  println(io, "\n", "═" ^ 88)
  println(io, @sprintf("  %s  |  %s", r.label, r.description))
  println(io, "  Ambient: ", r.ambient)
  println(io, "  HH dims: [", join(_fmt.(r.hh_dims), ", "), "]")
  println(io, "═" ^ 88)
  println(io)
  show(io, MIME("text/plain"), r.polyvector)
  println(io)
end

function main()
  println("\nComputing Hochschild cohomology of Küchle fourfolds...\n")

  results = HochschildKuechleResult[]
  t_total = 0.0

  for (label, k, n, weights, desc) in KUECHLE_FAMILIES
    print("  $label ... ")
    flush(stdout)
    t = @elapsed result = compute_hochschild_family(label, k, n, weights, desc)
    push!(results, result)
    t_total += t
    println(@sprintf("done  (%.2fs)", t))
  end

  println(@sprintf("\n  Total: %.2fs", t_total))
  println()
  show_summary(results)

  for r in results
    show_result(stdout, r)
  end
end

main()
