# Validate the exact case used for Mac/Linux reproducibility checks.
#
# Run from the package root with:
#   julia --project=. scripts/validate_exact_case.jl
#
# Expected current-code result:
#   score = 4.8783576387658645
#   chorus fraction = 0.320625
#   chorus count = 513
#   character count = 1087

using Printf

include(joinpath(@__DIR__, "..", "src", "SarfatiPelegChimera.jl"))

const EXPECTED = (
    score = 4.8783576387658645,
    synchrony_fraction = 0.320625,
    chorus_count = 513,
    char_count = 1087,
    freq_mean = 0.38068713011729466,
    freq_std = 0.0045542274862573765,
    freq_cv = 0.011963176913425362,
    lc_cx = 0.3950931921372677,
    lc_cy = 7.536484564072712,
    lc_points = 6000,
)

function check_close(name, got, expected; atol=1e-10)
    ok = isapprox(got, expected; atol=atol, rtol=0.0)
    @printf("%-22s got=% .16g expected=% .16g %s\n",
        name, got, expected, ok ? "OK" : "FAIL")
    return ok
end

function check_equal(name, got, expected)
    ok = got == expected
    @printf("%-22s got=%s expected=%s %s\n",
        name, string(got), string(expected), ok ? "OK" : "FAIL")
    return ok
end

function main()
    cfg = SarfatiPelegChimera.CFG
    case = (
        P = 9,
        kappa = 0.0485,
        tau0 = 1.185,
        sigtau = 0.062,
        ell = 12.0,
        delay_seed = 29,
        init_seed = 202,
    )

    println("Sarfati-Peleg exact-case validation")
    println("Julia version: ", VERSION)
    println("Threads: ", Threads.nthreads())
    println("Case: ", case)
    println()

    lc = SarfatiPelegChimera.build_limit_cycle(cfg)
    result = SarfatiPelegChimera.simulate_case(cfg, case, lc; keep_full=false)
    s = result.summary

    ok = true
    ok &= check_close("score", s.score, EXPECTED.score)
    ok &= check_close("chorus_fraction", s.synchrony_fraction, EXPECTED.synchrony_fraction)
    ok &= check_equal("chorus_count", s.chorus_count, EXPECTED.chorus_count)
    ok &= check_equal("char_count", s.char_count, EXPECTED.char_count)
    ok &= check_close("freq_mean", s.freq_mean, EXPECTED.freq_mean)
    ok &= check_close("freq_std", s.freq_std, EXPECTED.freq_std)
    ok &= check_close("freq_cv", s.freq_cv, EXPECTED.freq_cv)
    ok &= check_close("limit_cycle_cx", lc.cx, EXPECTED.lc_cx)
    ok &= check_close("limit_cycle_cy", lc.cy, EXPECTED.lc_cy)
    ok &= check_equal("limit_cycle_points", length(lc.x), EXPECTED.lc_points)

    println()
    if ok
        println("VALIDATION PASSED")
    else
        error("VALIDATION FAILED")
    end
end

main()
