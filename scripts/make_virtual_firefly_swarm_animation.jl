#!/usr/bin/env julia

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using DelimitedFiles
using Plots
using Random
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DATADIR = joinpath(ROOT, "results", "natural_chimera_case_01_gnuplot_data")
const OUTDIR = joinpath(ROOT, "figures", "virtual_firefly_swarm")
const N = 100
const TSTART = 260.0
const TSTOP = 330.0
const FPS = 3
const POSITION_SEED = 20261117
const LAYOUT_MODE = Symbol(get(ENV, "FIREFLY_LAYOUT", "random_field"))

mkpath(OUTDIR)

read_numeric_table(path) = Float64.(readdlm(path; comments=true, comment_char='#'))

# Recover the paper simulation and detect flashes exactly as in the diagnostics.
space = read_numeric_table(joinpath(DATADIR, "panel_01_space_time_y.dat"))
times = sort(unique(space[:, 1]))
nt = length(times)
Y = fill(NaN, nt, N)
time_index = Dict(t => k for (k, t) in enumerate(times))
for row in eachrow(space)
    Y[time_index[row[1]], Int(row[2])] = row[3]
end
late = findall(times .>= 260.0)
flash_threshold = quantile(vec(Y[late, :]), 0.82)
flash_times = [Float64[] for _ in 1:N]
for i in 1:N, k in 2:(nt - 1)
    if times[k] >= 260.0 && Y[k, i] > flash_threshold &&
       Y[k, i] > Y[k - 1, i] && Y[k, i] >= Y[k + 1, i]
        push!(flash_times[i], times[k])
    end
end

# Oscillator identity is encoded continuously by the measured character
# probability p_i: green (persistent chorus) -> gold (role switching) ->
# orange-red (persistent character).
probability_data = read_numeric_table(joinpath(DATADIR, "panel_11_character_probability.dat"))
character_probability = fill(0.5, N)
for row in eachrow(probability_data)
    character_probability[Int(row[1])] = clamp(row[2], 0.0, 1.0)
end
chorus_rgb = (0.30, 1.00, 0.38)
character_rgb = (1.00, 0.24, 0.045)
role_color(p) = RGB(
    (1 - p) * chorus_rgb[1] + p * character_rgb[1],
    (1 - p) * chorus_rgb[2] + p * character_rgb[2],
    (1 - p) * chorus_rgb[3] + p * character_rgb[3])
role_colors = role_color.(character_probability)
core_colors = [RGB(0.55 + 0.45c.r, 0.55 + 0.45c.g, 0.55 + 0.45c.b) for c in role_colors]

chorus_peak_data = read_numeric_table(joinpath(DATADIR, "panel_02_chorus_peaks.dat"))
chorus_times = chorus_peak_data[
    (chorus_peak_data[:, 1] .>= TSTART) .& (chorus_peak_data[:, 1] .<= TSTOP), 1]

# Reproducible virtual field: broad horizontal distribution, depth-dependent size,
# and small individual drift resembling freely moving fireflies.
Random.seed!(POSITION_SEED)
x0 = LAYOUT_MODE == :periodic_strip ? collect(1.0:N) : 4 .+ 92 .* rand(N)
y0 = 8 .+ 78 .* rand(N)
depth = 0.65 .+ 0.65 .* (1 .- y0 ./ 100) .+ 0.18 .* rand(N)
drift_ax = LAYOUT_MODE == :periodic_strip ? 0.03 .+ 0.07 .* rand(N) : 0.35 .+ 1.25 .* rand(N)
drift_ay = 0.20 .+ 0.75 .* rand(N)
drift_rate = 0.012 .+ 0.020 .* rand(N)
drift_phase_x = 2pi .* rand(N)
drift_phase_y = 2pi .* rand(N)

# Fixed faint stars and a layered grass silhouette give spatial depth without
# obscuring the firefly signal.
star_x = 100 .* rand(75)
star_y = 68 .+ 30 .* rand(75)
star_s = 0.35 .+ 0.9 .* rand(75)
grass_x = collect(0.0:0.35:100.0)
grass_y = 3.0 .+ 2.0 .* rand(length(grass_x)) .+ 2.0 .* sin.(0.18 .* grass_x).^2

function nearest_flash_distance(events, t)
    isempty(events) && return Inf
    return minimum(abs.(events .- t))
end

# A flash is rendered as a short, asymmetric pulse: rapid ignition and a slower
# fading glow. The 0.2 time-step data are sufficiently resolved for this effect.
function flash_brightness(events, t)
    isempty(events) && return 0.0
    previous = events[events .<= t]
    upcoming = events[events .> t]
    after = isempty(previous) ? Inf : t - previous[end]
    before = isempty(upcoming) ? Inf : upcoming[1] - t
    tail = after <= 1.15 ? exp(-after / 0.42) : 0.0
    rise = before <= 0.25 ? 0.32exp(-before / 0.10) : 0.0
    return clamp(max(tail, rise), 0.0, 1.0)
end

frame_indices = findall((times .>= TSTART) .& (times .<= TSTOP))

anim = @animate for k in frame_indices
    t = times[k]
    phase_t = t - TSTART
    x = x0 .+ drift_ax .* sin.(drift_rate .* phase_t .+ drift_phase_x)
    y = y0 .+ drift_ay .* sin.(0.83 .* drift_rate .* phase_t .+ drift_phase_y)
    b = [flash_brightness(flash_times[i], t) for i in 1:N]
    active = findall(b .> 0.035)
    strong = findall(b .> 0.50)

    p = plot(; xlims=(0, 100), ylims=(0, 100), size=(1000, 650),
        showaxis=LAYOUT_MODE == :periodic_strip ? :x : false,
        xticks=LAYOUT_MODE == :periodic_strip ? ([1, 25, 50, 75, 100], ["1", "25", "50", "75", "100"]) : :none,
        yticks=false,
        xlabel=LAYOUT_MODE == :periodic_strip ? "oscillator index i   (periodic boundary: 100 <-> 1)" : "",
        legend=false, background_color=RGB(0.008, 0.015, 0.055),
        foreground_color=:white, aspect_ratio=:equal,
        title=LAYOUT_MODE == :periodic_strip ?
            "SP-like firefly swarm on a periodic oscillator ring" :
            "A virtual SP-like firefly swarm",
        titlefontcolor=RGB(0.90, 0.93, 1.0), titlefontsize=14,
        margin=4Plots.mm)

    scatter!(p, star_x, star_y; markersize=star_s, markercolor=:white,
        markeralpha=0.28, markerstrokewidth=0, label=false)
    plot!(p, grass_x, grass_y; fillrange=0, fillcolor=RGB(0.01, 0.07, 0.035),
        fillalpha=1.0, linecolor=RGB(0.02, 0.12, 0.055), linewidth=1, label=false)

    # Dim bodies remain barely visible between flashes.
    scatter!(p, x, y; markersize=1.2 .* depth, markercolor=RGB(0.25, 0.34, 0.20),
        markeralpha=0.34, markerstrokewidth=0, label=false)

    if !isempty(active)
        # Three translucent layers produce a role-coloured bioluminescent halo.
        scatter!(p, x[active], y[active];
            markersize=(25 .+ 42 .* b[active]) .* depth[active],
            markercolor=role_colors[active], markeralpha=0.055,
            markerstrokewidth=0, label=false)
        scatter!(p, x[active], y[active];
            markersize=(12 .+ 24 .* b[active]) .* depth[active],
            markercolor=role_colors[active], markeralpha=0.14,
            markerstrokewidth=0, label=false)
        scatter!(p, x[active], y[active];
            markersize=(3 .+ 9 .* b[active]) .* depth[active],
            markercolor=role_colors[active], markeralpha=0.82,
            markerstrokecolor=core_colors[active], markerstrokewidth=0.25,
            label=false)
    end
    if !isempty(strong)
        scatter!(p, x[strong], y[strong];
            markersize=(1.5 .+ 3.2 .* b[strong]) .* depth[strong],
            markercolor=core_colors[strong], markeralpha=0.98,
            markerstrokewidth=0, label=false)
    end

    nbright = count(b .> 0.50)
    state_text = nbright >= 12 ? "collective chorus" : nbright > 0 ? "character flashes" : "between flashes"
    chorus_number = isempty(chorus_times) ? 0 : argmin(abs.(chorus_times .- t))
    at_chorus = chorus_number > 0 && abs(chorus_times[chorus_number] - t) <= 1.2
    time_color = at_chorus ? RGB(1.0, 0.93, 0.30) : RGB(0.72, 0.78, 0.90)
    time_size = at_chorus ? 21 : 18
    annotate!(p, 3.0, 96.0,
        text("t = $(round(t; digits=1))", time_size, time_color, :left))
    annotate!(p, 97.0, 96.0, text("chorus", 9, role_color(0.0), :right))
    annotate!(p, 97.0, 92.5, text("transition", 9, role_color(0.5), :right))
    annotate!(p, 97.0, 89.0, text("character", 9, role_color(1.0), :right))
    if at_chorus
        annotate!(p, 50.0, 91.5,
            text("CHORUS #$(chorus_number)", 18,
                 RGB(1.0, 0.90, 0.28), :center))
    end
    annotate!(p, 97.0, 5.5, text(state_text, 10,
        nbright >= 12 ? RGB(1.0, 0.88, 0.30) : RGB(0.62, 0.70, 0.82), :right))
    p
end

basename = LAYOUT_MODE == :periodic_strip ?
    "periodic_strip_sp_firefly_swarm" : "virtual_sp_firefly_swarm"
gif_path = joinpath(OUTDIR, basename * ".gif")
mp4_path = joinpath(OUTDIR, basename * ".mp4")
gif(anim, gif_path; fps=FPS)
mp4(anim, mp4_path; fps=FPS)

println("Frames: ", length(frame_indices), " at ", FPS, " fps")
println("Saved: ", gif_path)
println("Saved: ", mp4_path)
