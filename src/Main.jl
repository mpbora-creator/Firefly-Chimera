function main()
    println("Sarfati–Peleg-type chimera framework")
    println("Mode: ", RUN_MODE)
    println("Quick mode: ", QUICK)
    println("Julia threads available: ", nthreads())
    println("Delay history lookup: fractional linear interpolation")
    println("No A/B groups; random initial phases on limit cycle.")
    println("Sample manual case: ", SAMPLE_CASE)
    println(@sprintf("Score target: TARGET_CHORUS=%.2f, CHORUS_WIDTH=%.2f", TARGET_CHORUS, CHORUS_WIDTH))
    if RUN_MODE == :paper
        run_paper_mode()
        return
    end

    grid = active_grid()
    active_ells = active_ell_values(CFG, grid)
    println("Kernel: ", CFG.kernel)
    println("Active ell values: ", active_ells)
    println("Grid sizes: P=$(length(grid.P_values)), κ=$(length(grid.kappa_values)), τ0=$(length(grid.tau0_values)), στ=$(length(grid.sigtau_values)), ell=$(length(active_ells)), delay_seed=$(length(grid.delay_seeds)), init_seed=$(length(grid.init_seeds))")
    if RUN_MODE == :robustness
        c = robustness_center()
        println(@sprintf("Robustness centre: P=%d κ=%.6f τ0=%.6f στ=%.6f ell=%.2f", c.P, c.kappa, c.tau0, c.sigtau, c.ell))
        println("Robustness source best file: ", ROBUSTNESS_BEST_FILE, "  use_best_file=", ROBUSTNESS_USE_BEST_FILE)
        println(@sprintf("Robustness half-widths: ΔP=%d Δκ=%.6f Δτ0=%.6f Δστ=%.6f", ROBUST_DP, ROBUST_DKAPPA, ROBUST_DTAU0, ROBUST_DSIGTAU))
        println(@sprintf("Robustness success threshold: score >= %.3f", SUCCESS_SCORE_THRESHOLD))
    end

    println("Building single-firefly limit cycle once only; it is reused for every case...")
    lc = build_limit_cycle(CFG)
    println(@sprintf("Limit cycle ready: %d points, center=(%.6g, %.6g)", length(lc.x), lc.cx, lc.cy))

    cases = build_cases(CFG, grid)
    println("Total cases in grid: ", length(cases))

    ensure_live_csv(LIVE_CSV_PATH)
    logged_rows, completed_keys = RESUME_FROM_LIVE_CSV ? load_logged_rows(LIVE_CSV_PATH) : (Any[], Set{String}())
    if RESUME_FROM_LIVE_CSV
        println("Resume enabled. Previously completed cases in live CSV: ", length(completed_keys))
    end

    pending_indices = [idx for idx in eachindex(cases) if !(case_key(cases[idx]) in completed_keys)]
    println("Pending cases this run: ", length(pending_indices))
    if isempty(pending_indices)
        println("No pending cases. Regenerating summary from live CSV only.")
    end

    results = Vector{Any}(undef, length(pending_indices))
    progress_lock = ReentrantLock()
    log_lock = ReentrantLock()
    plot_lock = ReentrantLock()
    best_lock = ReentrantLock()
    completed_counter = Threads.Atomic{Int}(0)
    thread_counts = zeros(Int, Threads.maxthreadid())

    previous_best = isempty(logged_rows) ? -Inf : maximum([r.summary.score for r in logged_rows])
    best_score_ref = Ref(previous_best)
    if isfinite(previous_best)
        best_prev = sort(logged_rows, by = r -> r.summary.score, rev=true)[1].summary
        update_best_file(best_prev, BEST_PATH)
    end

    # Thread-parallel parameter scan with dynamic scheduling.
    # Each completed case is appended immediately to LIVE_CSV_PATH.
    # If the program crashes, restart it and completed cases are skipped.
    @threads :dynamic for outidx in eachindex(pending_indices)
        idx = pending_indices[outidx]
        case = cases[idx]
        tid = threadid()
        try
            res = simulate_case(CFG, case, lc; keep_full=false)
            results[outidx] = res
            s = res.summary

            # Count successful work on this thread.
            thread_counts[tid] += 1
            done = Threads.atomic_add!(completed_counter, 1) + 1

            # Crash-resilient append after every finished case.
            lock(log_lock) do
                append_live_csv(s, LIVE_CSV_PATH)
            end

            # Update best-so-far file immediately.
            lock(best_lock) do
                if s.score > best_score_ref[]
                    best_score_ref[] = s.score
                    update_best_file(s, BEST_PATH)
                end
            end

            # Save full diagnostics immediately for high-scoring cases.
            if SAVE_HIGH_SCORE_PLOTS && s.score >= HIGH_SCORE_PLOT_THRESHOLD
                full = simulate_case(CFG, case, lc; keep_full=true)
                outpng = high_score_plot_path(full.summary)
                lock(plot_lock) do
                    plot_case(full, 0, outpng)
                end
            end

            lock(progress_lock) do
                @printf("[%4d/%4d | done %4d/%4d] P=%2d κ=%.3f τ0=%.3f στ=%.3f ell=%.1f dseed=%d iseed=%d score=%7.3f chorus=%.2f char=%4d CV=%.3g thread=%d\n",
                    idx, length(cases), done, length(pending_indices), s.P, s.kappa, s.tau0, s.sigtau, s.ell,
                    s.delay_seed, s.init_seed, s.score, s.synchrony_fraction, s.char_count, s.freq_cv, tid)
                flush(stdout)
            end
        catch err
            thread_counts[tid] += 1
            done = Threads.atomic_add!(completed_counter, 1) + 1
            lock(progress_lock) do
                @warn "Case failed" idx case err tid done
            end
            results[outidx] = nothing
        end
    end

    println("Thread utilization (cases handled per Julia thread): ", thread_counts)
    println("Live CSV: ", LIVE_CSV_PATH)
    println("Best-so-far file: ", BEST_PATH)
    if SAVE_HIGH_SCORE_PLOTS
        println(@sprintf("High-score plots saved for score >= %.3f in %s", HIGH_SCORE_PLOT_THRESHOLD, HIGH_SCORE_DIR))
    end

    new_rows = [r for r in results if r !== nothing]
    rows = sort(vcat(logged_rows, new_rows), by = r -> r.summary.score, rev=true)
    csv_path = joinpath(OUTDIR, "$(OUT_PREFIX)_results.csv")
    summary_path = joinpath(OUTDIR, "$(OUT_PREFIX)_summary.png")
    write_csv(rows, csv_path)
    display(plot_summary(rows, summary_path))

    top20_path = joinpath(OUTDIR, "$(OUT_PREFIX)_top20.csv")
    write_top_cases(rows, top20_path; topn=20)
    println("Saved top-20 CSV: ", top20_path)

    if RUN_MODE == :robustness
        rob_rows = summarize_robustness(rows; success_threshold=SUCCESS_SCORE_THRESHOLD)
        rob_csv = joinpath(OUTDIR, "$(OUT_PREFIX)_robustness_summary.csv")
        write_robustness_summary(rob_rows, rob_csv)
        println("Saved robustness aggregate CSV: ", rob_csv)
        display(plot_robustness_heatmaps(rob_rows, OUTDIR; success_threshold=SUCCESS_SCORE_THRESHOLD))
        if !isempty(rob_rows)
            rb = rob_rows[1]
            @printf("Best robust aggregate: P=%d κ=%.5f τ0=%.5f στ=%.5f mean_score=%.3f success=%.2f n=%d\n",
                rb.P, rb.kappa, rb.tau0, rb.sigtau, rb.mean_score, rb.success_fraction, rb.n)
        end
    end

    println("\nTop cases with full diagnostic plots:")
    for (i,row) in enumerate(rows[1:min(CFG.top_cases_to_plot, length(rows))])
        full = simulate_case(CFG, row.case, lc; keep_full=true)
        out = joinpath(OUTDIR, OUT_PREFIX * @sprintf("_case_%02d.png", i))
        plot_case(full, i, out)

        #order_out = joinpath(OUTDIR, OUT_PREFIX * @sprintf("_case_%02d_order_verification.png", i))
        #plot_order_parameter_verification(full, i, order_out)

        #theta_out = joinpath(OUTDIR, OUT_PREFIX * @sprintf("_case_%02d_theta_tracks.png", i))

        #plot_theta_tracks_through_dip(
        #    full,
        #    i,
        #    theta_out;
        #    nshow = 16,
        #    half_window = 10.0
        #)

        s = full.summary
        @printf("case %02d -> %s | score=%.3f chorus=%.2f char=%d CV=%.3g\n",
            i, out, s.score, s.synchrony_fraction, s.char_count, s.freq_cv)
    end

    println("\nSaved CSV: ", csv_path)
    println("Saved summary: ", summary_path)
    println("Saved case figures: $(OUT_PREFIX)_case_XX.png")
end

