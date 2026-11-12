# -----------------------------
# Main
# -----------------------------

function robustness_groups(rows)
    groups = Dict{String,Vector{Any}}()
    for r in rows
        s = r.summary
        key = @sprintf("P=%d|k=%.8f|tau=%.8f|sig=%.8f|ell=%.8f", s.P, s.kappa, s.tau0, s.sigtau, s.ell)
        if !haskey(groups, key)
            groups[key] = Any[]
        end
        push!(groups[key], r)
    end
    return groups
end

function summarize_robustness(rows; success_threshold=SUCCESS_SCORE_THRESHOLD)
    out = NamedTuple[]
    for (_,rs) in robustness_groups(rows)
        ss = [r.summary for r in rs]
        scores = [s.score for s in ss]
        choruses = [s.synchrony_fraction for s in ss]
        cvs = [s.freq_cv for s in ss]
        dis = [s.di_track_score for s in ss]
        stables = [s.stable_character_fraction for s in ss]
        s0 = ss[1]
        push!(out, (;
            P=s0.P, kappa=s0.kappa, tau0=s0.tau0, sigtau=s0.sigtau, ell=s0.ell,
            n=length(ss),
            mean_score=mean(scores),
            std_score=length(scores) > 1 ? std(scores) : 0.0,
            max_score=maximum(scores),
            success_fraction=mean(scores .>= success_threshold),
            mean_chorus=mean(choruses),
            mean_freq_cv=mean(cvs),
            mean_di_track_score=mean(dis),
            mean_stable_character_fraction=mean(stables),
        ))
    end
    return sort(out, by = x -> x.mean_score, rev=true)
end

function write_robustness_summary(rows, path)
    cols = [:rank, :P, :kappa, :tau0, :sigtau, :ell, :n, :mean_score, :std_score, :max_score,
            :success_fraction, :mean_chorus, :mean_freq_cv, :mean_di_track_score,
            :mean_stable_character_fraction]
    open(path, "w") do io
        println(io, join(string.(cols), ","))
        for (rank,r) in enumerate(rows)
            vals = [rank, r.P, r.kappa, r.tau0, r.sigtau, r.ell, r.n, r.mean_score, r.std_score, r.max_score,
                    r.success_fraction, r.mean_chorus, r.mean_freq_cv, r.mean_di_track_score,
                    r.mean_stable_character_fraction]
            println(io, join(vals, ","))
        end
    end
end

function robustness_heatmap_matrix(summary_rows, Pfix, sigfix; field=:mean_score)
    kappas = sort(unique([r.kappa for r in summary_rows if r.P == Pfix && isapprox(r.sigtau, sigfix; atol=1e-12)]))
    taus = sort(unique([r.tau0 for r in summary_rows if r.P == Pfix && isapprox(r.sigtau, sigfix; atol=1e-12)]))
    Z = fill(NaN, length(taus), length(kappas))
    lookup = Dict{Tuple{Float64,Float64},Float64}()
    for r in summary_rows
        if r.P == Pfix && isapprox(r.sigtau, sigfix; atol=1e-12)
            lookup[(r.kappa, r.tau0)] = getfield(r, field)
        end
    end
    for (i,tau) in enumerate(taus), (j,k) in enumerate(kappas)
        if haskey(lookup, (k,tau))
            Z[i,j] = lookup[(k,tau)]
        end
    end
    return kappas, taus, Z
end

function plot_robustness_heatmaps(summary_rows, outdir; success_threshold=SUCCESS_SCORE_THRESHOLD)
    isempty(summary_rows) && return nothing
    mkpath(outdir)
    best = summary_rows[1]
    Pfix = best.P
    sigfix = best.sigtau

    kappas, taus, Zscore = robustness_heatmap_matrix(summary_rows, Pfix, sigfix; field=:mean_score)
    _, _, Zsucc = robustness_heatmap_matrix(summary_rows, Pfix, sigfix; field=:success_fraction)

    if isempty(kappas) || isempty(taus)
        return nothing
    end

    p1 = heatmap(kappas, taus, Zscore; xlabel="κ", ylabel="τ0",
        title=@sprintf("Mean score, P=%d, στ=%.4f", Pfix, sigfix),
        colorbar_title="mean score")
    p2 = heatmap(kappas, taus, Zsucc; xlabel="κ", ylabel="τ0",
        title=@sprintf("Success fraction, score ≥ %.2f", success_threshold),
        colorbar_title="success")

    fig = plot(p1, p2; layout=(1,2), size=(1500,650))
    path = joinpath(outdir, @sprintf("robustness_heatmaps_P%d_sigtau_%.4f.png", Pfix, sigfix))
    savefig(fig, path)
    return fig
end
