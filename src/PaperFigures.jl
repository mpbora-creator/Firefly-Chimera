function write_top_cases(rows, path; topn=20)
    top = rows[1:min(topn, length(rows))]
    write_csv(top, path)
end

function read_rows_from_csv(path::AbstractString)
    rows, _ = load_logged_rows(path)
    return sort(rows, by = r -> r.summary.score, rev=true)
end

function run_paper_mode()
    println("Paper mode: regenerate diagnostics from top completed search cases.")
    println("Source CSV: ", PAPER_SOURCE_CSV)
    rows = read_rows_from_csv(PAPER_SOURCE_CSV)
    if isempty(rows)
        println("No parseable completed rows found in paper source CSV.")
        return
    end
    mkpath(PAPER_OUTDIR)
    n = min(PAPER_TOP_N, length(rows))
    println("Top cases available: ", length(rows), "; plotting top ", n)
    table_path = joinpath(PAPER_OUTDIR, "paper_top$(n)_table.csv")
    write_top_cases(rows, table_path; topn=n)
    println("Saved paper top-case table: ", table_path)

    println("Building single-firefly limit cycle once only for paper diagnostics...")
    lc = build_limit_cycle(CFG)
    println(@sprintf("Limit cycle ready: %d points, center=(%.6g, %.6g)", length(lc.x), lc.cx, lc.cy))

    for (i,row) in enumerate(rows[1:n])
        full = simulate_case(CFG, row.case, lc; keep_full=true)
        out = joinpath(PAPER_OUTDIR, @sprintf("paper_case_%02d_score_%.3f_P%d_k%.4f_tau%.4f_sig%.4f.png",
            i, full.summary.score, full.summary.P, full.summary.kappa, full.summary.tau0, full.summary.sigtau))
        plot_case(full, i, out)
        s = full.summary
        @printf("paper case %02d -> %s | score=%.3f chorus=%.2f stable=%d d_i=%.3f CV=%.3g\n",
            i, out, s.score, s.synchrony_fraction, s.stable_character_count, s.di_track_score, s.freq_cv)
    end
    println("Paper figures saved in: ", PAPER_OUTDIR)
end
