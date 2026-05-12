# ═══════════════════════════════════════════════════════════════════════════════
#  CuriousGrassmanian.jl — Hochschild Cohomology for OGr(n-1,2n+1), n ≥ 4
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties

function main()
  for n in 4:8
    a = n - 1
    b=2 * n + 1
    X = OGr(a, b)
    T = filtered_tangent_bundle(X)
    flag = true
    for q in 1:rank_bundle(T)
      if has_higher_cohomology(exterior_power(T, q))
        msg = "Hochschild cohomology is not global for OGr($a, $b) " *
              "in degree q = $q"
        println(msg)
        flag = false
        break
      end
    end
    if flag
      msg =
        "Not enough information to conclude non-global Hochschild " *
        "cohomology for OGr($a, $b)"
      println(msg)
    end
  end
end

main()
