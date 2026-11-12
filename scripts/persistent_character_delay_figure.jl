# Auxiliary figure derived exclusively from the original numerical model and
# original Diagnostics.jl definitions. No analytical/validation module is loaded.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

module OriginalSarfatiPelegNumerics
include(joinpath(@__DIR__, "..", "src", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "Utilities.jl"))
include(joinpath(@__DIR__, "..", "src", "LimitCycle.jl"))
include(joinpath(@__DIR__, "..", "src", "Coupling.jl"))
include(joinpath(@__DIR__, "..", "src", "Simulation.jl"))
include(joinpath(@__DIR__, "..", "src", "Diagnostics.jl"))
end

using .OriginalSarfatiPelegNumerics
using Statistics, DelimitedFiles, Plots
const S = OriginalSarfatiPelegNumerics

case = (P=9, kappa=0.0485, tau0=1.185, sigtau=0.062, ell=12.0,
        delay_seed=29, init_seed=202)
lc = S.build_limit_cycle(S.CFG)
result = S.simulate_case(S.CFG, case, lc; keep_full=true)
tr = result.delay_tracks
ids = tr.stable_oscillators
means = tr.stable_delay_mean
absmedian = median(abs.(means))

outdir = joinpath(S.FIGURES_DIR, "persistent_character_delay")
mkpath(outdir)

# Save the exact plotted quantities for reproducibility.
open(joinpath(outdir, "stable_character_tracks.csv"), "w") do io
    println(io, "oscillator,stable_run_mean,stable_run_std,run_length,start_cycle,stop_cycle")
    for q in eachindex(ids)
        println(io, join((ids[q], means[q], tr.stable_delay_std[q],
                          tr.stable_run_lengths[q], tr.stable_start_cycles[q],
                          tr.stable_stop_cycles[q]), ','))
    end
end

p1 = plot(xlabel="chorus cycle n", ylabel="dᵢ⁽ⁿ⁾",
          title="Persistent character-delay tracks", legend=false)
for (q, i) in enumerate(ids)
    rr = tr.stable_start_cycles[q]:tr.stable_stop_cycles[q]
    vals = tr.delay_matrix[i, rr]
    plot!(p1, tr.cycles[rr], vals; color=:gray, alpha=0.34, linewidth=1)
end
highlight = findall(x -> isapprox(x, absmedian; atol=1e-10), means)
for q in highlight
    i = ids[q]
    rr = tr.stable_start_cycles[q]:tr.stable_stop_cycles[q]
    plot!(p1, tr.cycles[rr], tr.delay_matrix[i, rr]; color=:crimson,
          linewidth=2.5, marker=:circle, markersize=2.5)
end
hline!(p1, [absmedian]; color=:black, linestyle=:dash, linewidth=2,
       label=false)
annotate!(p1, (1.4, absmedian + 0.07,
              text("median |μᵢ| = $(round(absmedian, digits=7))", 9, :left)))

p2 = scatter(ids, means; xlabel="oscillator i", ylabel="stable-run mean μᵢ",
             title="Signed persistent-character means", legend=false,
             markercolor=:steelblue, markersize=4)
scatter!(p2, ids[highlight], means[highlight]; markercolor=:crimson,
         markersize=6, markerstrokewidth=0)
hline!(p2, [absmedian]; color=:black, linestyle=:dash, linewidth=2)
annotate!(p2, (55, absmedian + 0.055, text("0.5016667", 9, :left)))

bins = collect(0.28:0.05:0.88)
p3 = histogram(means; bins=bins, xlabel="stable-run mean μᵢ",
               ylabel="number of tracks", title="Distribution of signed means",
               legend=false, color=:steelblue, alpha=0.75)
vline!(p3, [absmedian]; color=:black, linestyle=:dash, linewidth=2)

caption = "Original Case 1 numerics: 30 persistent tracks; all selected signed means are positive in this realization."
layout = @layout [a{0.48h}; b c]
fig = plot(p1, p2, p3; layout=layout, size=(1100, 900), margin=7Plots.mm,
           plot_title=caption, plot_titlefontsize=11)
path = joinpath(outdir, "persistent_character_delay_evidence.png")
savefig(fig, path)
println("Figure: ", path)
println("Stable tracks: ", length(ids))
println("Signed range: ", extrema(means))
println("Median absolute stable-run mean: ", absmedian)
