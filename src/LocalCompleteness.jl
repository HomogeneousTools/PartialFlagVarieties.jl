export is_locally_complete

"""
    is_locally_complete(Z::ZeroLocus) -> Union{Bool,Nothing}

Check a sufficient criterion for the family of sections of the defining bundle
to be locally complete at a smooth regular zero locus `Z`.

Returns `true` when both cohomology groups
``\\mathrm{H}^1(Z, \\mathrm{T}_X|_Z)`` and
``\\mathrm{H}^1(X, \\mathcal{E} \\otimes \\mathcal{I}_Z)`` vanish, where
``X`` is the ambient variety and ``\\mathcal{E}`` is its defining bundle.
Otherwise returns `nothing`: the sufficient criterion is inconclusive,
so the family is not asserted to be incomplete.

The restriction and normal sequences give the Kodaira–Spencer map

```math
\\begin{equation}
\\mathrm{H}^0(X, \\mathcal{E})
\\longrightarrow \\mathrm{H}^0(Z, \\mathcal{E}|_Z)
\\longrightarrow \\mathrm{H}^1(Z, \\mathrm{T}_Z).
\\end{equation}
```

The two vanishings force both arrows to be surjective, so the
[Kodaira–Spencer completeness theorem](https://doi.org/10.1007/s13373-015-0070-1)
(Corollary 157) applies to the smooth base of regular sections.
The ideal-sheaf twist is computed from the Koszul resolution with its
``\\mathcal{E}`` term omitted, in the full ambient cohomological range.
The actual filtered tangent bundle is used on a general partial flag variety.

No chosen section, generic maximal-rank assumption, or unobstructedness
assumption is required for this certificate.
It applies to every smooth regular section for which these vanishings are
forced by the recorded data.
Smoothness and regularity of the zero locus are assumed, as for
[`tangent_bundle(::ZeroLocus)`](@ref).

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> Z = zero_locus(line_bundle(projective_space(4), 5));

julia> is_locally_complete(Z)
true

julia> K3 = zero_locus(line_bundle(projective_space(3), 4));

julia> is_locally_complete(K3) === nothing
true
```
"""
function is_locally_complete(Z::ZeroLocus)
  X = ambient_variety(Z)
  counter = Ref(0)
  if dimension(Z) > 0
    tangent_h1 = _restrict_to_zero_locus_les(Z, filtered_tangent_bundle(X), counter)[2]
    is_zero_expr(tangent_h1) || return nothing
  end

  codimension(Z) == 0 && return true
  # E ⊗ ∧^r E* → ⋯ → E ⊗ E* resolves E ⊗ I_Z, a sheaf on X.
  # Do not truncate its cohomology to dim(Z), even when Z is a point.
  terms = _koszul_dimensions(Z, defining_bundle(Z))[2:end]
  ideal_h1 = solve_koszul_filtration_symbolic(terms, dimension(X), counter)[1]
  is_zero_expr(ideal_h1) ? true : nothing
end
