# ═══════════════════════════════════════════════════════════════════════════════
#  PartialFlagVarieties.jl
#
#  A Julia package for computing with partial flag varieties G/P using
#  compile-time specialization via Lie.jl.
#
#  Provides:
#   - MarkedDynkinType: type-level encoding of G/P
#   - PartialFlagVariety: named user-facing wrapper
#   - IrrepLevi: irreducible representations of the Levi subgroup
#   - CompletelyReducibleBundle: equivariant bundles on G/P
#   - Cohomology: sheaf cohomology via Borel–Weil–Bott
#   - Constructions: named constructors (Gr, OGr, SGr, ℙⁿ, etc.)
# ═══════════════════════════════════════════════════════════════════════════════

module PartialFlagVarieties

using Lie
using Preferences
using StaticArrays
using LinearAlgebra
using Combinatorics

# ─── Extend Lie.jl functions (avoids name collisions) ────────────────────────

import Lie: dimension, dual, tensor_product, exterior_power, symmetric_power
import Lie: n_components

# ─── Core types and infrastructure ───────────────────────────────────────────

include("Lie.jl")
include("MarkedDynkinType.jl")
include("PartialFlagVariety.jl")
include("IrrepLevi.jl")
include("CompletelyReducibleBundle.jl")
include("FilteredBundle.jl")
include("UniversalBundles.jl")
include("Cohomology.jl")
include("Hochschild.jl")
include("Constructions.jl")
include("Koszul.jl")
include("ZeroLoci.jl")

# ─── Reexport commonly used Lie.jl types ─────────────────────────────────────

using Lie: TypeA, TypeB, TypeC, TypeD, TypeE, TypeF4, TypeG2
using Lie: DynkinType, SimpleDynkinType, ProductDynkinType
using Lie: WeightLatticeElem, WeylCharacter, fundamental_weight
using Lie: rank, degree

export TypeA, TypeB, TypeC, TypeD, TypeE, TypeF4, TypeG2
export DynkinType, SimpleDynkinType, ProductDynkinType
export WeightLatticeElem, WeylCharacter, fundamental_weight

# ─── Export Lie-extended functions and package-specific functions ─────────────

export dimension, dual, tensor_product, exterior_power, symmetric_power
export n_components
export tangent_weights, positive_nonparabolic_roots, positive_parabolic_roots
export direct_sum, twist
export cartan_type, cartan_type_with_ordering, parse_dynkin_type
export graded_pieces, total_bundle, filtered_tangent_bundle, n_filtration_steps
export universal_subbundle, universal_quotient_bundle, spinor_bundle
export tautological_bundles, quotient_bundles
export hodge_numbers, twisted_hodge_numbers, hochschild_cohomology
export PolyvectorParallelogram
export ZeroLocus, zero_locus, ambient_variety, defining_bundle
export codimension, normal_bundle, conormal_bundle
export koszul_terms, cohomology_on_restriction
export is_calabi_yau, is_calabi_yau_candidate
export solve_ses_cohomology, solve_koszul_filtration

# ─── Startup banner ──────────────────────────────────────────────────────────

function _print_banner()
  v = pkgversion(@__MODULE__)
  version_str = v === nothing ? "dev" : string(v)

  println()
  printstyled(" ██████╗     "; color=:blue)
  print("██╗")
  printstyled("██████╗  "; color=:red)
  println(" │  equivariant bundles on partial flag varieties G/P")

  printstyled("██╔════╝    "; color=:blue)
  print("██╔╝")
  printstyled("██╔══██╗"; color=:red)
  println("  │  cohomology via the Borel–Weil–Bott theorem")

  printstyled("██║  ███╗  "; color=:blue)
  print("██╔╝ ")
  printstyled("██████╔╝"; color=:red)
  println("  │")

  printstyled("██║   ██║ "; color=:blue)
  print("██╔╝  ")
  printstyled("██╔═══╝ "; color=:red)
  println("  │  Docs:    https://homogeneous.tools")

  printstyled("╚██████╔╝"; color=:blue)
  print("██╔╝   ")
  printstyled("██║     "; color=:red)
  println("  │  Version: ", version_str)

  printstyled(" ╚═════╝ "; color=:blue)
  print("╚═╝    ")
  printstyled("╚═╝     "; color=:red)
  println("  │")
end

# ─── Initialization ──────────────────────────────────────────────────────────

function __init__()
  # Suppress the Lie.jl startup banner: PartialFlagVarieties will show its own
  set_preferences!(Lie, "show_banner" => false)

  if displaysize(stdout)[2] >= 80
    _print_banner()
  end
  return nothing
end

if Base.VERSION >= v"1.9"
  # Type A
  let X = Gr(2, 4)
    precompile(dimension, (typeof(X),))
    precompile(betti_numbers, (typeof(X),))
  end
  precompile(levi_type, (Type{MarkedDynkinType{TypeA{3},(2,)}},))
  precompile(decomposition_matrix, (Type{MarkedDynkinType{TypeA{3},(2,)}},))

  # Type B
  let X = partial_flag_variety(TypeB{3}, (1,))
    precompile(dimension, (typeof(X),))
  end

  # Type D
  let X = partial_flag_variety(TypeD{5}, (5,))
    precompile(dimension, (typeof(X),))
  end
end

end # module PartialFlagVarieties
