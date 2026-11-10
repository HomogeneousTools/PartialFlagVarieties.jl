export local_completeness

"""
    local_completeness(Z::ZeroLocus) -> NamedTuple

Check a sufficient criterion for the family of sections of the defining bundle
to be locally complete at a smooth regular zero locus `Z`.

Returns `status = :complete` when both cohomology groups
``\\mathrm{H}^1(Z, \\mathrm{T}_X|_Z)`` and
``\\mathrm{H}^1(X, \\mathcal{E} \\otimes \\mathcal{I}_Z)`` vanish, where
``X`` is the ambient variety and ``\\mathcal{E}`` is its defining bundle.
Otherwise returns `status = :undetermined`, which does not assert that the
family is incomplete.
The fields `ambient_tangent_h1` and `ideal_twist_h1` contain the corresponding
dimension-valued [`AffineExpr`](@ref) entries, possibly symbolic.

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

julia> local_completeness(Z)
(status = :complete, ambient_tangent_h1 = 0, ideal_twist_h1 = 0)

julia> K3 = zero_locus(line_bundle(projective_space(3), 4));

julia> local_completeness(K3)
(status = :undetermined, ambient_tangent_h1 = 1, ideal_twist_h1 = 0)
```
"""
function local_completeness(Z::ZeroLocus)
  X = ambient_variety(Z)
  counter = Ref(0)
  tangent_h1 = if dimension(Z) == 0
    AffineExpr(0)
  else
    _restrict_to_zero_locus_les(Z, filtered_tangent_bundle(X), counter)[2]
  end

  ideal_h1 = if codimension(Z) == 0
    AffineExpr(0)
  else
    # E ⊗ ∧^r E* → ⋯ → E ⊗ E* resolves E ⊗ I_Z, a sheaf on X.
    # Do not truncate its cohomology to dim(Z), even when Z is a point.
    terms = _koszul_dimensions(Z, defining_bundle(Z))[2:end]
    solve_koszul_filtration_symbolic(terms, dimension(X), counter)[1]
  end

  entries = [tangent_h1, ideal_h1]
  _renumber_variables!(entries)
  (
    status=all(is_zero_expr, entries) ? :complete : :undetermined,
    ambient_tangent_h1=entries[1],
    ideal_twist_h1=entries[2],
  )
end
