# ═══════════════════════════════════════════════════════════════════════════════
#  BottVanishingE6.jl — Bott vanishing on maximal parabolic of E6
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties

function main()
    for i in 1:6
        DT = MarkedDynkinType(TypeE{6}, (i,))
        X = PartialFlagVariety(DT)
        Ω = dual(filtered_tangent_bundle(X))
        flag = true
        for q in 1:rank_bundle(Ω)
            for n in 1:20
                Et = tensor_product(
                    exterior_power(Ω, q),
                    line_bundle(Ω.variety, n)
                )
                if has_higher_cohomology(Et)
                    msg = "Bott vanishing fails for exterior power q = $q " *
                          "on $X with twist n = $n"
                    println(msg)
                    flag = false
                    break
                end
            end
            if !flag
                break
            end
        end
        if flag
            println("Not enough evidence to disprove Bott vanishing for $X")
        end
    end
end

main() 