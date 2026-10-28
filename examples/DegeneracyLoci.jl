using PartialFlagVarieties

# Determinantal fourfolds with product and relative Grassmann presentations.
function determinantal_examples()
  X = projective_space(8)
  OX, L = structure_sheaf(X), line_bundle(X, 1)
  [
    degeneracy_locus(4 * OX, 4 * L, 2),
    degeneracy_locus(direct_sum(2 * OX, dual(L)), 3 * L, 1; presentation=:relative),
  ]
end

# All five rows of Benedetti--Filippini--Manivel--Tanturri, Table 8.
# The odd symplectic IGr(2,7) is a hyperplane section of Gr(2,7).
function orbital_fano_fourfolds()
  X = Gr(3, 6)
  first_row = partially_decomposable_locus(
    direct_sum(dual(universal_subbundle(X)), 3 * structure_sheaf(X))
  )
  Y = Gr(2, 7)
  second_row = partially_decomposable_locus(
    direct_sum(universal_quotient_bundle(Y), structure_sheaf(Y));
    base_equations=line_bundle(Y, 1),
  )
  third_row = partially_decomposable_locus(
    direct_sum(dual(universal_subbundle(Y)), 4 * structure_sheaf(Y));
    base_equations=line_bundle(Y, 1),
  )

  # The projectivized rank-one nilpotent orbit is Fl(1,n-1;n), embedded by
  # O(1,1) in P(sl_n). The Q¹³ row cuts it by one quadric; the P²⁰ row
  # cuts Fl(1,4;5) by three hyperplanes in P²³.
  A, B = flag_variety(4, [1, 3]), flag_variety(5, [1, 4])
  fourth_row = zero_locus(line_bundle(A, [2, 2]))
  fifth_row = zero_locus(3 * line_bundle(B, [1, 1]))
  [first_row, second_row, third_row, fourth_row, fifth_row]
end

function table8_euler_invariants(Z)
  Ω = cotangent_bundle(Z)
  (;
    degree=degree(anticanonical_bundle(Z)),
    chi1=chi(Ω),
    chi2=chi(exterior_power(Ω, 2)),
    anticanonical_sections=cohomology(anticanonical_bundle(Z))[0],
  )
end

if abspath(PROGRAM_FILE) == @__FILE__
  configure_caches!(; budget=512 * 1024^2)
  for (i, Z) in enumerate(orbital_fano_fourfolds())
    println("Table 8, row ", i, ": ", table8_euler_invariants(Z))
  end
end
