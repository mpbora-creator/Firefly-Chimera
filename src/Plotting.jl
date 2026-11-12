function write_csv(rows, path)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(CSV_HEADER), ","))
        for (r,row) in enumerate(rows)
            println(io, join(csv_values_from_summary(row.summary; rank=r), ","))
        end
    end
end


# ===== GNUPLOT DATA EXPORT HELPERS =====
# These routines dump the numerical data behind every diagnostic panel into
# plain text files. Each file starts with commented Gnuplot commands.
# They do not alter any plot or diagnostic; they only write auxiliary data.

function _gpstr(x)
    if x === missing
        return "NaN"
    elseif x isa Number
        if isfinite(float(x))
            return @sprintf("%.12g", float(x))
        else
            return "NaN"
        end
    else
        return string(x)
    end
end

function _gp_write_header(io, title, command_lines)
    println(io, "# ", title)
    println(io, "#")
    println(io, "# Suggested Gnuplot command(s):")
    for line in command_lines
        println(io, "# ", line)
    end
    println(io, "#")
end

function _gp_write_xy(path, title, colnames, rows, command_lines)
    mkpath(dirname(path))
    open(path, "w") do io
        _gp_write_header(io, title, command_lines)
        println(io, "# columns: ", join(colnames, " "))
        for row in rows
            println(io, join(_gpstr.(row), " "))
        end
    end
end

function _gp_write_text(path, title, lines)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", title)
        println(io, "# This is a text-only diagnostics panel.")
        println(io, "# Suggested Gnuplot usage: use set label commands manually, or inspect as text.")
        println(io, "#")
        for line in lines
            println(io, "# ", line)
        end
    end
end

function _gp_zvalue(Z, iy, ix, ny, nx)
    if size(Z,1) == ny && size(Z,2) == nx
        return Z[iy, ix]
    elseif size(Z,1) == nx && size(Z,2) == ny
        return Z[ix, iy]
    else
        return NaN
    end
end

function _gp_write_matrix(path, title, xvals, yvals, Z, command_lines; xname="x", yname="y", zname="z")
    mkpath(dirname(path))
    nx = length(xvals)
    ny = length(yvals)
    open(path, "w") do io
        _gp_write_header(io, title, command_lines)
        println(io, "# columns: ", xname, " ", yname, " ", zname)
        for iy in 1:ny
            for ix in 1:nx
                z = _gp_zvalue(Z, iy, ix, ny, nx)
                println(io, join((_gpstr(xvals[ix]), _gpstr(yvals[iy]), _gpstr(z)), " "))
            end
            println(io)
        end
    end
end

function _gp_dump_case_data(res, idx, figpath, lt, late_cols, text_lines)
    base = splitext(figpath)[1] * "_gnuplot_data"
    mkpath(base)

    # Panel 01
    _gp_write_matrix(joinpath(base, @sprintf("panel_%02d_space_time_y.dat", 1)),
        @sprintf("Panel 01: Case %02d y_i(t) space-time", idx),
        res.times, 1:CFG.N, res.Ys,
        ["set view map",
         "set xlabel 'time'; set ylabel 'oscillator index'",
         "set title 'Case y_i(t) space-time'",
         "splot 'panel_01_space_time_y.dat' using 1:2:3 with pm3d notitle"];
        xname="time", yname="oscillator_index", zname="y")

    # Panel 02
    rows2 = [(res.chorus_bins[k], res.chorus_counts[k]) for k in eachindex(res.chorus_bins)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_flash_count.dat", 2)),
        "Panel 02: collective flash count", ["time", "flashes_per_bin"], rows2,
        ["set xlabel 'time'; set ylabel 'flashes/bin'",
         "plot 'panel_02_flash_count.dat' using 1:2 with lines lw 2 title 'flash count', \\",
         "     'panel_02_chorus_peaks.dat' using 1:2 with points pt 7 ps 1.2 title 'chorus peaks'"])
    chorus_y = [res.chorus_counts[argmin(abs.(res.chorus_bins .- tc))] for tc in res.chorus]
    rows2b = [(res.chorus[k], chorus_y[k]) for k in eachindex(res.chorus)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_chorus_peaks.dat", 2)),
        "Panel 02 overlay: detected chorus peaks", ["time", "peak_flash_count"], rows2b,
        ["plot 'panel_02_chorus_peaks.dat' using 1:2 with points pt 7 ps 1.2 title 'chorus peaks'"])

    # Panel 03
    rows3 = [(lt[k], res.R1[k], res.R2[k], mean(res.R1), mean(res.R2)) for k in eachindex(lt)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_global_order.dat", 3)),
        "Panel 03: global order parameters", ["time", "R1", "R2", "mean_R1", "mean_R2"], rows3,
        ["set xlabel 'time'; set ylabel 'global order'",
         "plot 'panel_03_global_order.dat' using 1:2 with lines title 'R1', \\",
         "     '' using 1:3 with lines title 'R2', \\",
         "     '' using 1:4 with lines lw 2 title '<R1>', \\",
         "     '' using 1:5 with lines lw 2 title '<R2>'"])

    # Panel 04
    _gp_write_matrix(joinpath(base, @sprintf("panel_%02d_local_order_Z.dat", 4)),
        "Panel 04: local order Z_i(t)", lt, 1:CFG.N, res.Zmat,
        ["set view map",
         "set xlabel 'time'; set ylabel 'oscillator index'",
         "splot 'panel_04_local_order_Z.dat' using 1:2:3 with pm3d notitle"];
        xname="time", yname="oscillator_index", zname="Z")

    # Panel 05
    rows5 = [(i, res.final_theta[i]) for i in 1:CFG.N]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_final_phase.dat", 5)),
        "Panel 05: final reconstructed phase", ["oscillator_index", "theta"], rows5,
        ["set xlabel 'oscillator index'; set ylabel 'wrapped phase'",
         "plot 'panel_05_final_phase.dat' using 1:2 with points pt 7 title 'theta_i'"])

    # Panel 06
    rows6 = [(i, res.final_Z[i]) for i in 1:CFG.N]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_final_local_coherence.dat", 6)),
        "Panel 06: final local coherence", ["oscillator_index", "Z_i"], rows6,
        ["set xlabel 'oscillator index'; set ylabel 'Z_i'",
         "plot 'panel_06_final_local_coherence.dat' using 1:2 with points pt 7 title 'Z_i'"])

    # Panel 07
    rows7 = [(x,) for x in res.classif.character_delays if isfinite(x)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_character_delays_raw.dat", 7)),
        "Panel 07: character delays raw values", ["delay"], rows7,
        ["binwidth=0.1",
         "bin(x,w)=w*floor(x/w)+w/2.0",
         "set boxwidth binwidth",
         "plot 'panel_07_character_delays_raw.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title 'character delays'"])

    # Panel 08
    rows8 = [(x,) for x in res.omega if isfinite(x)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_frequency_raw.dat", 8)),
        "Panel 08: mean phase frequency raw values", ["omega"], rows8,
        ["binwidth=0.0005",
         "bin(x,w)=w*floor(x/w)+w/2.0",
         "set boxwidth binwidth",
         "plot 'panel_08_frequency_raw.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title 'frequency profile'"])

    # Panel 09
    _gp_write_matrix(joinpath(base, @sprintf("panel_%02d_phase_difference.dat", 9)),
        "Panel 09: phase difference Delta_i(t)", lt, 1:CFG.N, res.Dmat,
        ["set view map",
         "set xlabel 'time'; set ylabel 'oscillator index'",
         "splot 'panel_09_phase_difference.dat' using 1:2:3 with pm3d notitle"];
        xname="time", yname="oscillator_index", zname="Delta")

    # Panel 10
    _gp_write_matrix(joinpath(base, @sprintf("panel_%02d_cycle_delay.dat", 10)),
        "Panel 10: cycle-resolved delay d_i(n)", res.delay_tracks.cycles, 1:CFG.N, res.delay_tracks.delay_matrix,
        ["set view map",
         "set xlabel 'chorus cycle n'; set ylabel 'oscillator index'",
         "splot 'panel_10_cycle_delay.dat' using 1:2:3 with pm3d notitle"];
        xname="cycle_n", yname="oscillator_index", zname="d_i")

    # Panel 11
    rows11 = [(i, res.delay_tracks.character_probability[i]) for i in 1:CFG.N]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_character_probability.dat", 11)),
        "Panel 11: character-role probability p_i", ["oscillator_index", "p_i"], rows11,
        ["set xlabel 'oscillator index'; set ylabel 'p_i'",
         "plot 'panel_11_character_probability.dat' using 1:2 with points pt 7 title 'p_i', \\",
         "     0.2 with lines dashtype 2 notitle, 0.8 with lines dashtype 2 notitle"])

    # Panel 12
    rows12 = [(x,) for x in res.delay_tracks.character_probability if isfinite(x)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_character_probability_raw.dat", 12)),
        "Panel 12: character-role probability raw values", ["p_i"], rows12,
        ["binwidth=0.05",
         "bin(x,w)=w*floor(x/w)+w/2.0",
         "set boxwidth binwidth",
         "plot 'panel_12_character_probability_raw.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title 'p_i distribution'"])

    # Panel 13
    _gp_write_matrix(joinpath(base, @sprintf("panel_%02d_role_raster.dat", 13)),
        "Panel 13: role raster", res.delay_tracks.cycles, 1:CFG.N, res.delay_tracks.character_state_matrix,
        ["set view map",
         "set xlabel 'chorus cycle n'; set ylabel 'oscillator index'",
         "splot 'panel_13_role_raster.dat' using 1:2:3 with pm3d notitle"];
        xname="cycle_n", yname="oscillator_index", zname="role_state")

    # Panel 14
    rows14 = [("chorus", res.summary.persistent_chorus_count),
              ("switcher", res.summary.switching_oscillator_count),
              ("character", res.summary.persistent_character_count)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_persistent_role_counts.dat", 14)),
        "Panel 14: persistent role counts", ["role", "count"], rows14,
        ["set style data histograms",
         "set style fill solid",
         "set xlabel 'oscillator role'; set ylabel 'count'",
         "plot 'panel_14_persistent_role_counts.dat' using 2:xtic(1) title 'count'"])

    # Panel 15
    _gp_write_text(joinpath(base, @sprintf("panel_%02d_diagnostics_text.dat", 15)),
        "Panel 15: diagnostics summary", text_lines)

    # Panel 16
    rows16 = [(i, res.delay_tracks.role_switch_count[i]) for i in 1:CFG.N]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_role_switch_count.dat", 16)),
        "Panel 16: role-switch count", ["oscillator_index", "switch_count"], rows16,
        ["set xlabel 'oscillator index'; set ylabel 'switch count'",
         "plot 'panel_16_role_switch_count.dat' using 1:2 with points pt 7 title 'switch count'"])

    # Panel 17
    rows17 = [(i, res.delay_variability[i]) for i in 1:CFG.N]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_delay_variability.dat", 17)),
        "Panel 17: delay variability sigma_i", ["oscillator_index", "sigma_i"], rows17,
        ["set xlabel 'oscillator index'; set ylabel 'sigma_i'",
         "plot 'panel_17_delay_variability.dat' using 1:2 with points pt 7 title 'sigma_i'"])

    # Panel 18
    _gp_write_matrix(joinpath(base, @sprintf("panel_%02d_delay_residual.dat", 18)),
        "Panel 18: delay residual d_i(n)-<d_i>", res.delay_tracks.cycles, 1:CFG.N, res.delay_residual_matrix,
        ["set view map",
         "set xlabel 'chorus cycle n'; set ylabel 'oscillator index'",
         "splot 'panel_18_delay_residual.dat' using 1:2:3 with pm3d notitle"];
        xname="cycle_n", yname="oscillator_index", zname="delay_residual")

    # Panel 19
    _gp_write_matrix(joinpath(base, @sprintf("panel_%02d_delay_covariance_matrix.dat", 19)),
        "Panel 19: delay covariance C_ij", 1:CFG.N, 1:CFG.N, res.delay_covariance_matrix,
        ["set view map",
         "set xlabel 'oscillator index j'; set ylabel 'oscillator index i'",
         "splot 'panel_19_delay_covariance_matrix.dat' using 1:2:3 with pm3d notitle"];
        xname="oscillator_j", yname="oscillator_i", zname="C_ij")

    # Panel 20
    rows20 = [(res.delay_cov_distance.distances[k], res.delay_cov_distance.mean_cov[k]) for k in eachindex(res.delay_cov_distance.distances)]
    _gp_write_xy(joinpath(base, @sprintf("panel_%02d_delay_covariance_vs_distance.dat", 20)),
        "Panel 20: delay covariance versus ring distance", ["ring_distance", "mean_C_r"], rows20,
        ["P=" * string(res.case.P),
         "set xlabel 'ring distance |i-j|'; set ylabel 'mean C(r)'",
         "plot 'panel_20_delay_covariance_vs_distance.dat' using 1:2 with linespoints title 'C(r)', \\",
         "     P with lines dashtype 2 title 'P'"])

    return base
end
# ===== END GNUPLOT DATA EXPORT HELPERS =====


function plot_summary(rows, path)
    mkpath(dirname(path))
    scores = [r.summary.score for r in rows]
    p1 = histogram(scores; bins=25, xlabel="score", ylabel="count", title="Sarfati–Peleg similarity score", legend=false)
    p2 = scatter([r.summary.synchrony_fraction for r in rows], [r.summary.score for r in rows];
        xlabel="chorus fraction", ylabel="score", title=@sprintf("Chorus target %.2f ± %.2f", TARGET_CHORUS, CHORUS_WIDTH), legend=false)
    p3 = scatter([r.summary.freq_cv for r in rows], [r.summary.score for r in rows];
        xlabel="frequency CV", ylabel="score", title="Common-frequency locking", legend=false)
    p4 = plot(legend=false, grid=false, axis=false, ticks=false, title="Top candidates")
    lines = String[]
    for (i,r) in enumerate(rows[1:min(12,length(rows))])
        s = r.summary
        push!(lines, @sprintf("%2d. P=%d κ=%.3f τ0=%.3f στ=%.3f seed=(%d,%d) score=%.3f chorus=%.2f Cs=%.2f d_i=%.2f CV=%.3g",
            i, s.P, s.kappa, s.tau0, s.sigtau, s.delay_seed, s.init_seed, s.score,
            s.synchrony_fraction, s.chorus_score, s.di_track_score, s.freq_cv))
    end
    annotate!(p4, 0.02, 0.98, text(join(lines,"\n"), 8, :left, :top))
    fig = plot(p1,p2,p3,p4; layout=(2,2), size=(1500,900))
    savefig(fig, path)
    return fig
end

function plot_case(res, idx, path)
    mkpath(dirname(path))
    times = res.times
    late_cols = findall(times .>= CFG.t_transient)
    lt = times[late_cols]

    p1 = heatmap(times, 1:CFG.N, res.Ys; xlabel="time", ylabel="oscillator index",
        title=@sprintf("Case %02d: y_i(t) space-time", idx), colorbar_title="y")

    p2 = plot(res.chorus_bins, res.chorus_counts; xlabel="time", ylabel="flashes/bin",
        title="Collective flash count and detected chorus", label="flash count", lw=2)
    if !isempty(res.chorus)
    chorus_y = [
        res.chorus_counts[argmin(abs.(res.chorus_bins .- tc))]
        for tc in res.chorus
    ]

    scatter!(
        p2,
        res.chorus,
        chorus_y;
        label = "chorus peaks",
        markercolor = :black,
        markersize = 4
    )
end

    # Global order panel: controlled by hardcoded GLOBAL_ORDER_PLOT_MODE in Diagnostics.jl.
    # This changes only the display, not res.R1/res.R2 or any score.
    R1_mean_line = fill(mean(res.R1), length(lt))
    R2_mean_line = fill(mean(res.R2), length(lt))

    if GLOBAL_ORDER_PLOT_MODE == 0
        p3 = plot(lt, res.R1; label="R1", xlabel="time", ylabel="global order",
            title="Global order parameters", lw=2)
        plot!(p3, lt, res.R2; label="R2", lw=2)
    elseif GLOBAL_ORDER_PLOT_MODE == 1
        p3 = plot(lt, R1_mean_line; label="<R1>", xlabel="time", ylabel="mean global order",
            title="Mean global order parameters", lw=2)
        plot!(p3, lt, R2_mean_line; label="<R2>", lw=2)
        annotate!(p3, lt[1], 0.05, text(@sprintf("<R1>=%.3f\n<R2>=%.3f", mean(res.R1), mean(res.R2)), 8, :left, :bottom))
    elseif GLOBAL_ORDER_PLOT_MODE == 2
        p3 = plot(lt, res.R1; label="R1", xlabel="time", ylabel="global order",
            title="Global order parameters with means", lw=1, alpha=0.35)
        plot!(p3, lt, res.R2; label="R2", lw=1, alpha=0.35)
        plot!(p3, lt, R1_mean_line; label="<R1>", lw=3)
        plot!(p3, lt, R2_mean_line; label="<R2>", lw=3)
    else
        error("Invalid GLOBAL_ORDER_PLOT_MODE=$(GLOBAL_ORDER_PLOT_MODE). Use 0, 1, or 2.")
    end

    p4 = heatmap(lt, 1:CFG.N, res.Zmat; xlabel="time", ylabel="oscillator index",
        title="Local order Z_i(t)", colorbar_title="Z")

    p5 = scatter(1:CFG.N, res.final_theta; xlabel="oscillator index", ylabel="wrapped phase",
        title="Final reconstructed phase", legend=false, markersize=4)

    p6 = scatter(1:CFG.N, res.final_Z; xlabel="oscillator index", ylabel="Z_i",
        title="Final local coherence", legend=false, markersize=4)

    p7 = histogram(res.classif.character_delays; bins=30, xlabel="delay from nearest chorus", ylabel="count",
        title="Character delays", legend=false)

    p8 = histogram(res.omega; bins=24, xlabel="mean phase frequency", ylabel="count",
        title=@sprintf("Frequency profile, CV=%.3g", res.summary.freq_cv), legend=false)

    p9 = heatmap(lt, 1:CFG.N, res.Dmat; xlabel="time", ylabel="oscillator index",
        title="Phase difference Δ_i(t)", colorbar_title="Δ")

    delay_plot_matrix = nan_to_missing_matrix(res.delay_tracks.delay_matrix)
    p10 = heatmap(res.delay_tracks.cycles, 1:CFG.N, delay_plot_matrix;
        xlabel="chorus cycle n", ylabel="oscillator index",
        title="Cycle-resolved delay d_i(n)", colorbar_title="d_i")

    pchar = res.delay_tracks.character_probability
    valid_pchar = [x for x in pchar if isfinite(x)]
    p11 = scatter(1:CFG.N, pchar; xlabel="oscillator index", ylabel="p_i",
        title="Character-role probability p_i", legend=false, markersize=4, ylim=(-0.05,1.05))
    hline!(p11, [0.2, 0.8]; linestyle=:dash, label=false)

    p12 = histogram(valid_pchar; bins=20, xlabel="p_i", ylabel="count",
        title="Distribution of character-role probability", legend=false, xlim=(0,1))

    p13 = heatmap(res.delay_tracks.cycles, 1:CFG.N, res.delay_tracks.character_state_matrix;
        xlabel="chorus cycle n", ylabel="oscillator index", title="Role raster: chorus ↔ character",
        colorbar_title="0 chorus, 1 char", clims=(0,1))

    role_names = ["chorus", "switcher", "character"]
    role_counts = [res.summary.persistent_chorus_count, res.summary.switching_oscillator_count, res.summary.persistent_character_count]
    p14 = bar(role_names, role_counts; xlabel="oscillator role", ylabel="count",
        title="Persistent role counts", legend=false)

    # ===== ROLE-SWITCH DYNAMICS DIAGNOSTIC PLOT =====
    p16 = scatter(1:CFG.N, res.delay_tracks.role_switch_count;
        xlabel="oscillator index", ylabel="switch count",
        title="Role-switch count", legend=false, markersize=4)
    # ===== END ROLE-SWITCH DYNAMICS DIAGNOSTIC PLOT =====

    # ===== DELAY-PERSISTENCE DIAGNOSTIC PLOTS =====
    # 1. Delay variability: sigma_i = std_n[d_i(n)].
    delay_var = res.delay_variability
    p17 = scatter(1:CFG.N, delay_var;
        xlabel="oscillator index", ylabel="σ_i",
        title="Delay variability σ_i", legend=false, markersize=4)

    # 2. Delay residual map: d_i(n) - <d_i>.
    delay_resid_plot_matrix = nan_to_missing_matrix(res.delay_residual_matrix)
    p18 = heatmap(res.delay_tracks.cycles, 1:CFG.N, delay_resid_plot_matrix;
        xlabel="chorus cycle n", ylabel="oscillator index",
        title="Delay residual d_i(n)-⟨d_i⟩", colorbar_title="residual")

    # 3. Delay covariance matrix: C_ij = cov_n(d_i,d_j).
    delay_cov_plot_matrix = nan_to_missing_matrix(res.delay_covariance_matrix)
    p19 = heatmap(1:CFG.N, 1:CFG.N, delay_cov_plot_matrix;
        xlabel="oscillator index j", ylabel="oscillator index i",
        title="Delay covariance C_ij", colorbar_title="cov")
    # ===== END DELAY-PERSISTENCE DIAGNOSTIC PLOTS =====

    stable_run_text = hasproperty(res.summary, :stable_character_mean_run) ? res.summary.stable_character_mean_run : NaN

    text_lines = [
        @sprintf("Sarfati-Peleg similarity score = %.4f", res.summary.score),
        @sprintf("P=%d, κ=%.4f, τ0=%.4f, στ=%.4f", res.summary.P, res.summary.kappa, res.summary.tau0, res.summary.sigtau),
        @sprintf("target chorus = %.2f ± %.2f", res.summary.target_chorus, res.summary.chorus_width),
        @sprintf("chorus fraction = %.3f, chorus score = %.3f", res.summary.synchrony_fraction, res.summary.chorus_score),
        @sprintf("freq score = %.3f, d_i score = %.3f", res.summary.common_freq_score, res.summary.di_track_score),
        @sprintf("chorus flashes = %d, character flashes = %d", res.summary.chorus_count, res.summary.char_count),
        @sprintf("character delay median/span = %.4f / %.4f", res.summary.delay_median, res.summary.delay_span),
        @sprintf("stable delayed oscillators = %d", res.summary.stable_character_count),
        @sprintf("mean |d_i| / std(d_i) = %.4f / %.4f", res.summary.stable_character_mean_delay, res.summary.stable_character_mean_std),
        @sprintf("mean stable run length = %.2f cycles", stable_run_text),
        @sprintf("p_i mean/std = %.3f / %.3f", res.summary.mean_character_probability, res.summary.std_character_probability),
        @sprintf("role counts: chorus=%d, switcher=%d, character=%d", res.summary.persistent_chorus_count, res.summary.switching_oscillator_count, res.summary.persistent_character_count),
        @sprintf("role switches: mean=%.2f, max=%d, mean lifetime=%.2f cycles", res.summary.mean_role_switch_count, res.summary.max_role_switch_count, res.summary.mean_role_lifetime),
        @sprintf("R1=%.3f, R2=%.3f", res.summary.R1_mean, res.summary.R2_mean),
        @sprintf("Z mean/std = %.3f / %.3f", res.summary.Z_mean, res.summary.Z_std),
        @sprintf("delay σ mean/std = %.3f / %.3f", res.summary.delay_variability_mean, res.summary.delay_variability_std),
        @sprintf("delay cov mean/std = %.3f / %.3f", res.summary.delay_covariance_mean, res.summary.delay_covariance_std),
        @sprintf("C(r): C0=%.3f, C1=%.3f, C(P)=%.3f", res.summary.delay_covariance_distance0, res.summary.delay_covariance_distance1, res.summary.delay_covariance_distanceP),
        @sprintf("freq mean/std/CV = %.4f / %.4f / %.4g", res.summary.freq_mean, res.summary.freq_std, res.summary.freq_cv),
        "Interpretation: p_i≈0 persistent chorus; p_i≈1 persistent character;",
        "intermediate p_i indicates switching between roles."
    ]
    p15 = plot(legend=false, grid=false, axis=false, ticks=false, title="Diagnostics")
    annotate!(p15, 0.02, 0.98, text(join(text_lines,"\n"), 7, :left, :top))

    # 4. Mean covariance versus ring distance: C(r)=<C_ij>_{r_ij=r}.
    # r=0 is the diagonal mean variance; r>0 shows spatial correlation
    # of delay fluctuations across the ring.
    p20 = plot(res.delay_cov_distance.distances, res.delay_cov_distance.mean_cov;
        xlabel="ring distance |i-j|", ylabel="mean C(r)",
        title="Delay covariance vs distance", legend=false, lw=2, marker=:circle)
    vline!(p20, [res.case.P]; linestyle=:dash, label=false)

    # Dump all panel data for external Gnuplot reproduction.
    _gp_dump_case_data(res, idx, path, lt, late_cols, text_lines)

    fig = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20; layout=(5,4), size=(2400,2350))
    savefig(fig, path)
    return fig
end
