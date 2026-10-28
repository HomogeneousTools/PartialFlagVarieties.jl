# A dominant graded piece is globally generated and has no higher cohomology.
# Induction through the filtration therefore lifts its generating sections.
_generation_certificate(F::CompletelyReducibleBundle) = is_globally_generated(F)
_generation_certificate(F::FilteredBundle) = all(is_globally_generated, graded_pieces(F))

function _require_generation(F::_AmbientBundle, name::AbstractString)
  _generation_certificate(F) && return nothing
  reason = if F isa FilteredBundle
    "Global generation could not be certified"
  else
    "The bundle is not globally generated"
  end
  throw(ArgumentError("$reason for $name."))
end

function _check_degeneracy_geometry(
  section::_AmbientBundle, cuts::_AmbientBundle, lower_codimension::Integer;
  check_global_generation::Bool=false, check_smoothness::Bool=false,
)
  variety(cuts) == variety(section) || throw(ArgumentError("Incorrect base equations."))
  if check_global_generation || check_smoothness
    _require_generation(section, "the section bundle")
    _require_generation(cuts, "the base equations")
  end
  if check_smoothness
    d = dimension(variety(section)) - rank(cuts)
    lower_codimension > d || throw(
      ArgumentError("A deeper stratum may occur; smoothness could not be certified.")
    )
  end
  nothing
end

function _relative_generation_certificate(F::RelativeGrassmannianBundle)
  G = variety(F)
  b, n = rank(dynkin_type(base_variety(G))), Int(rank(G.bundle))
  for piece in _relative_pieces(F.representation)
    # The central GL(n) character is unrestricted; the other fibre coordinates
    # must be dominant for the relative evaluation map to be surjective.
    for rep in components(piece)
      c = coefficients(p_dominant_weight(rep))
      all(>=(0), c[(b + 1):(b + n - 1)]) || return false
    end
    images = pushforward(G, RelativeGrassmannianBundle(G, piece))
    _generation_certificate(images[0]) || return false
  end
  true
end

function _relative_ample_certificate(L::RelativeGrassmannianBundle)
  rank(L) == 1 || return false
  G = variety(L)
  if G.k == rank(G.bundle)
    return is_ample_line_bundle(_relative_total(pushforward(G, L)[0]))
  end
  R = _relative_total(L.representation)
  c = coefficients(p_dominant_weight(only(components(R))))
  c[rank(dynkin_type(base_variety(G))) + G.k] > 0 || return false
  H = line_bundle(base_variety(G), ones(Int, picard_rank(base_variety(G))))
  # Relative very ampleness plus generation gives a closed embedding into
  # X × P^N. Tensoring with an ample line on X is then ample.
  _relative_generation_certificate(tensor_product(L, pullback(G, dual(H))))
end

"""
    is_strongly_fano(Z::RelativeZeroLocus)

A sufficient certificate that the adjunction anticanonical line has an ample
lift to the relative ambient. A false result means the certificate failed;
it does not prove that the locus is not Fano. Regularity and nonemptiness
remain assumptions about the section.
"""
function is_strongly_fano(Z::RelativeZeroLocus)
  G = ambient_variety(Z)
  L = tensor_product(anticanonical_bundle(G), dual(det(defining_bundle(Z))))
  _relative_ample_certificate(L)
end

function _positive_polarization(G::GrassmannBundle)
  X = base_variety(G)
  H = line_bundle(X, ones(Int, picard_rank(X)))
  G.k == rank(G.bundle) && return pullback(G, H)
  images = pushforward(G, O(G, 1))
  # Twist every graded direct-image summand into the dominant chamber.
  t = 0
  for piece in _relative_pieces(images[0])
    for rep in components(piece), i in marked_nodes(X)
      t = max(t, -coefficients(p_dominant_weight(rep))[i])
    end
  end
  tensor_product(O(G, 1), pullback(G, _line_power(H, max(0, t) + 1)))
end

function _check_degeneracy_nonempty(
  Z::Union{ZeroLocus,RelativeZeroLocus}; check_nonempty::Bool=false
)
  check_nonempty || return Z
  A = ambient_variety(Z)
  L = if A isa GrassmannBundle
    _positive_polarization(A)
  else
    line_bundle(A, ones(Int, picard_rank(A)))
  end
  volume = degree(restrict(Z, L))
  volume > 0 || throw(
    ArgumentError(
      if iszero(volume)
        "The expected regular locus is empty."
      else
        "A negative polarization degree contradicts the regularity assumptions."
      end,
    ),
  )
  Z
end
