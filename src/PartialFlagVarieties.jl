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
using StaticArrays
using LinearAlgebra
using Combinatorics

# ─── Extend Lie.jl functions (avoids name collisions) ────────────────────────

import Lie: dimension, dual, tensor_product, exterior_power, symmetric_power
import Lie: n_components

# ─── Core types and infrastructure ───────────────────────────────────────────

include("MarkedDynkinType.jl")
include("PartialFlagVariety.jl")
include("IrrepLevi.jl")
include("CompletelyReducibleBundle.jl")
include("Cohomology.jl")
include("Constructions.jl")

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

# ─── Precompilation hints ────────────────────────────────────────────────────

# Warm up key @generated methods for common types
function __init__()
  # Intentionally left minimal; Julia 1.9+ pkgimage precompilation
  # handles @generated functions automatically.
end

if Base.VERSION >= v"1.9"
  # Type A
  precompile(dimension, (Type{MarkedDynkinType{TypeA{3},(2,)}},))
  precompile(levi_type, (Type{MarkedDynkinType{TypeA{3},(2,)}},))
  precompile(betti_numbers, (Type{MarkedDynkinType{TypeA{3},(2,)}},))
  precompile(special_matrix, (Type{MarkedDynkinType{TypeA{3},(2,)}},))

  # Type B
  precompile(dimension, (Type{MarkedDynkinType{TypeB{3},(1,)}},))

  # Type D
  precompile(dimension, (Type{MarkedDynkinType{TypeD{5},(5,)}},))
end

end # module PartialFlagVarieties
