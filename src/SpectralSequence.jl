# ═══════════════════════════════════════════════════════════════════════════════
#  SpectralSequence.jl — Spectral Sequences
# ═══════════════════════════════════════════════════════════════════════════════

export is_spectral_sequence_degenerate
# ═══════════════════════════════════════════════════════════════════════════════
#  Helper functions
# ═══════════════════════════════════════════════════════════════════════════════

"""
    check_degeneracy_spectral_sequence(E::Matrix{Int}) -> Bool

Check if the spectral sequence represented by the E₁ page will degenerate at
some page. Returns false if there exists a possible differential between
non-zero entries, true otherwise.
"""
function check_degeneracy_spectral_sequence(E::Matrix{Int})
  d, n = size(E)
  for i in 1:(d-1), j in 1:(n-1)
    if E[i, j] != 0 && any(!=(0), @view E[i+1, j+1:n])
      return false
    end
  end
  return true
end
# ═══════════════════════════════════════════════════════════════════════════════
#  Spectral Sequence
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _construct_E1_page(F::FilteredBundle) -> Matrix{WeylCharacter}

Construct the E₁ page of the spectral sequence from a filtered bundle.
Helper function shared by `is_spectral_sequence_degenerate` and `E1_page`.
"""
function _construct_E1_page(F::FilteredBundle)
  n = n_filtration_steps(F)
  bundles = reverse(graded_pieces(F)) # Reverse to match spectral sequence convention
  d = dimension(variety(F))
  DT = dynkin_type(variety(F))
  R = rank(variety(F))
  E = Matrix{WeylCharacter{DT,R}}(undef, d, n)

  for j in 1:n
    H = cohomology(bundles[j])
    for i in 1:d
      E[i, j] = H[i-1]
    end
  end
  return E
end

function is_spectral_sequence_degenerate(F::FilteredBundle)
  E = _construct_E1_page(F)
  d, n = size(E)
  
  # Extract all weights that appear in the spectral sequence
  weights = Set{WeightLatticeElem}()
  for i in 1:d, j in 1:n
    if E[i, j] != 0
      union!(weights, keys(E[i, j].terms))
    end
  end

  # Check degeneracy for each weight's isotropic component
  for weight in weights
    E_iso = zeros(Int, d, n)
    for i in 1:d, j in 1:n
      if (char = E[i, j]) != 0 && haskey(char.terms, weight)
        E_iso[i, j] = char.terms[weight]
      end
    end
    if !check_degeneracy_spectral_sequence(E_iso)
      return false
    end
  end
  return true
end