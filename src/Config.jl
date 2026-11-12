# ============================================================
# Natural Sarfati–Peleg-type firefly chimera search
# ------------------------------------------------------------
# Goal:
#   Search for spontaneously emerging mixed chorus/character states
#   in N identical two-variable relaxation firefly oscillators.
#
# Key design choices:
#   * No imposed A/B groups.
#   * Random initial phases on the single-firefly limit cycle.
#   * Identical intrinsic oscillator parameters for all fireflies.
#   * Homogeneous nonlocal distance-weighted delayed coupling.
#   * Phase is reconstructed only after integration.
#   * Chorus/characters are detected from simulated flashes, following
#     the experimental logic of Sarfati & Peleg.
#
# Run search:
#   julia -t 10 sarfati_peleg_chimera_modular.jl
#
# Run robustness/sensitivity mode:
#   CHIMERA_MODE=robustness julia -t 10 sarfati_peleg_chimera_modular.jl
#
# Quick test:
#   CHIMERA_QUICK=1 julia -t 10 sarfati_peleg_chimera_modular.jl
#
# Outputs:
#   natural_chimera_results.csv
#   natural_chimera_summary.png
#   natural_chimera_case_01.png, ...
# ============================================================

using Random
using Base.Threads
using Statistics
using LinearAlgebra
using Printf
using DelimitedFiles
using Plots

const QUICK = get(ENV, "CHIMERA_QUICK", "0") == "1"

# One code base, three drivers:
#   CHIMERA_MODE=search      -> ordinary parameter search
#   CHIMERA_MODE=robustness  -> local robustness/sensitivity study around a chosen/best case
#   CHIMERA_MODE=paper       -> regenerate diagnostic plots and table from top cases only
const RUN_MODE = Symbol(get(ENV, "CHIMERA_MODE", "search"))
if !(RUN_MODE in (:search, :robustness, :paper))
    error("CHIMERA_MODE must be search, robustness, or paper, got $(RUN_MODE)")
end

# Directory layout
# Config.jl lives in src/, so the project root is one level above @__DIR__.
const SRCDIR = @__DIR__
const ROOTDIR = dirname(SRCDIR)
const RESULTS_DIR = joinpath(ROOTDIR, "results")
const FIGURES_DIR = joinpath(ROOTDIR, "figures")
mkpath(RESULTS_DIR)
mkpath(FIGURES_DIR)

# Most tabular/search outputs go to results/. Diagnostic images go to figures/.
const OUTDIR = RESULTS_DIR
const OUT_PREFIX = RUN_MODE == :robustness ? "robustness_chimera" :
                   RUN_MODE == :paper ? "paper_chimera" : "natural_chimera"
const PROGRESS_EVERY = QUICK ? 1 : 50

# -----------------------------
# Crash-resilient logging / resume
# -----------------------------
# The live CSV is appended after every completed case and flushed immediately.
# If the run is interrupted, restart the same program; completed cases are skipped.
const LIVE_CSV_PATH = joinpath(RESULTS_DIR, "$(OUT_PREFIX)_results_live.csv")
const BEST_PATH     = joinpath(RESULTS_DIR, "$(OUT_PREFIX)_best_so_far.txt")
const HIGH_SCORE_DIR = joinpath(FIGURES_DIR, "$(OUT_PREFIX)_high_score_cases")
# Resume is opt-in for reproducibility runs.  Older live CSV files may have
# been produced by an earlier diagnostic definition, so the default is to
# recompute rather than silently reuse stale classifications.
const RESUME_FROM_LIVE_CSV = get(ENV, "CHIMERA_RESUME", "0") == "1"
const SAVE_HIGH_SCORE_PLOTS = false
const HIGH_SCORE_PLOT_THRESHOLD = 4.0

# -----------------------------
# Sarfati–Peleg similarity-score target
# -----------------------------
# Tune these two first if you want the search to prefer a larger or smaller chorus.
# TARGET_CHORUS is the preferred fraction of flashes classified as chorus.
# CHORUS_WIDTH controls how tolerant the score is around TARGET_CHORUS.
const TARGET_CHORUS = 0.65
const CHORUS_WIDTH  = 0.15

# Relative weights in the new Sarfati–Peleg similarity score.
# Larger values make the corresponding feature more important.
const W_CHORUS      = 1.60
const W_FREQUENCY   = 1.35
const W_DELAY_SPAN  = 1.00
const W_DELAY_OFFSET= 0.75
const W_DI_TRACKS   = 1.25
const W_LOCAL_MIX   = 0.90
const W_CHAR_SYNC   = 0.45
const W_ACTIVITY    = 0.35

# -----------------------------
# Cycle tracking for Sarfati-Peleg characters
# -----------------------------
# These constants control the d_i(n) analysis.  A persistent "character"
# is an oscillator that remains off-chorus with nearly constant delay for
# several consecutive chorus cycles.  This is closer to the Sarfati-Peleg
# trajectory criterion than counting isolated off-beat flashes.
const CHARACTER_TRACK_MIN_RUN = 5          # minimum consecutive/off-chorus cycles
const CHARACTER_TRACK_STD_MAX = 0.22       # maximum std of d_i(n) within a stable run
const CHARACTER_TRACK_MIN_DELAY = 0.24     # exclude near-chorus oscillators
const CHARACTER_TRACK_MAX_GAP = 1          # allow this many missed chorus cycles inside a run

# -----------------------------
# User-tunable search parameters
# -----------------------------
Base.@kwdef struct SearchConfig
    # Network
    N::Int = 100
    geometry::Symbol = :ring             # :ring or :disk
    kernel::Symbol = :ring_uniform       # :ring_uniform or :exponential

    # Intrinsic oscillator: dx/dt = x^2 - beta*x*y; dy/dt = eta*x + rho*(x-1)*y
    beta::Float64 = 0.052
    rho::Float64 = 1.8
    eta::Float64 = 8.2

    # Integration
    dt::Float64 = QUICK ? 0.02 : 0.01
    t_end::Float64 = QUICK ? 180.0 : 520.0
    t_transient::Float64 = QUICK ? 80.0 : 260.0
    sample_dt::Float64 = QUICK ? 0.20 : 0.20

    # Limit cycle construction
    lc_relax_time::Float64 = 250.0
    lc_dt::Float64 = 0.01
    lc_keep_points::Int = 6000

    # Flash detection
    flash_variable::Symbol = :y          # :x or :y
    flash_quantile::Float64 = 0.82       # threshold from late-time signal quantile
    chorus_bin_dt::Float64 = 0.25        # bins for collective flash count
    chorus_window::Float64 = 0.18        # flashes within this time of chorus are chorus
    min_chorus_peak_fraction::Float64 = 0.12

    # Diagnostics
    local_order_delta::Int = 6
    top_cases_to_plot::Int = QUICK ? 2 : 8
    rng_seed::Int = 12345
end

const CFG = SearchConfig()

# -----------------------------
# Manual search region
# -----------------------------
# Edit these arrays freely. The first row below is deliberately narrow and fast.
# Use broader arrays once the program is working on your machine.
const SEARCH_GRID = (; 
    P_values       = QUICK ? [10, 16, 20] : [9],       # ring half-neighbours; P/N ~ 0.12-0.24
    kappa_values   = QUICK ? [0.05, 0.07, 0.09] : collect(0.048:0.0001:0.049),
    tau0_values    = QUICK ? [1.1, 1.2] : collect(1.175:0.001:1.19),
    sigtau_values  = QUICK ? [0.01, 0.03, 0.05] : collect(0.05:0.001:0.07),
    # ell is used only for kernel = :exponential.  For the default
    # kernel = :ring_uniform, keep ell_values as a singleton to avoid
    # duplicate runs that differ only by an unused ell value.
    ell_values     = [12.0],
    delay_seeds    = QUICK ? [11] : [11, 17, 29],
    init_seeds     = QUICK ? [101] : [101, 202]
)

# A single recommended starting case for manual exploration.
const SAMPLE_CASE = (; P=16, kappa=0.070, tau0=1.10, sigtau=0.04, ell=12.0, delay_seed=11, init_seed=101)

# -----------------------------
# Local robustness/sensitivity region
# -----------------------------
# Robustness mode can either use this manual centre, or automatically read
# the latest search best-so-far file.  The automatic file is normally
# natural_chimera_best_so_far.txt, generated by search mode.
#
# Activate robustness by:
#   CHIMERA_MODE=robustness julia -t 10 sarfati_peleg_chimera_modular.jl
#
# Activate paper plotting by:
#   CHIMERA_MODE=paper julia -t 10 sarfati_peleg_chimera_modular.jl

const ROBUSTNESS_USE_BEST_FILE = get(ENV, "ROBUST_FROM_BEST", "1") == "1"
const ROBUSTNESS_BEST_FILE = joinpath(RESULTS_DIR, get(ENV, "ROBUST_BEST_FILE", "natural_chimera_best_so_far.txt"))

# Manual fallback centre. Edit these if you do not want auto-centering.
const ROBUST_P       = parse(Int,     get(ENV, "ROBUST_P", "9"))
const ROBUST_KAPPA   = parse(Float64, get(ENV, "ROBUST_KAPPA", "0.0485"))
const ROBUST_TAU0    = parse(Float64, get(ENV, "ROBUST_TAU0", "1.185"))
const ROBUST_SIGTAU  = parse(Float64, get(ENV, "ROBUST_SIGTAU", "0.062"))
const ROBUST_ELL     = parse(Float64, get(ENV, "ROBUST_ELL", "12.0"))

# Local search widths around the robustness centre.
const ROBUST_DP       = parse(Int,     get(ENV, "ROBUST_DP", "1"))
const ROBUST_DKAPPA   = parse(Float64, get(ENV, "ROBUST_DKAPPA", "0.0015"))
const ROBUST_DTAU0    = parse(Float64, get(ENV, "ROBUST_DTAU0", "0.005"))
const ROBUST_DSIGTAU  = parse(Float64, get(ENV, "ROBUST_DSIGTAU", "0.005"))

const ROBUST_STEP_KAPPA  = parse(Float64, get(ENV, "ROBUST_STEP_KAPPA", "0.0005"))
const ROBUST_STEP_TAU0   = parse(Float64, get(ENV, "ROBUST_STEP_TAU0", "0.0025"))
const ROBUST_STEP_SIGTAU = parse(Float64, get(ENV, "ROBUST_STEP_SIGTAU", "0.0025"))

const ROBUST_DELAY_SEEDS = QUICK ? [29] : [11, 17, 29, 37, 43]
const ROBUST_INIT_SEEDS  = QUICK ? [202] : [101, 202, 303]

# Paper mode settings. The paper driver reads top cases from the search live CSV
# unless PAPER_SOURCE_CSV is specified. It does not run a new search grid.
const PAPER_SOURCE_CSV = joinpath(RESULTS_DIR, get(ENV, "PAPER_SOURCE_CSV", "natural_chimera_results_live.csv"))
const PAPER_TOP_N = parse(Int, get(ENV, "PAPER_TOP_N", "20"))
const PAPER_OUTDIR = joinpath(FIGURES_DIR, get(ENV, "PAPER_OUTDIR", "paper_figures"))

const SUCCESS_SCORE_THRESHOLD = 5.5

function parse_best_so_far(path::AbstractString)
    vals = Dict{Symbol,Any}()
    if !isfile(path)
        return nothing
    end
    for line in readlines(path)
        if occursin("=", line)
            parts = split(line, "="; limit=2)
            key = Symbol(strip(parts[1]))
            val = strip(parts[2])
            try
                if key in (:P, :delay_seed, :init_seed, :char_count)
                    vals[key] = parse_int_like(val)
                else
                    vals[key] = parse(Float64, val)
                end
            catch
                # ignore non-numeric lines
            end
        end
    end
    needed = (:P, :kappa, :tau0, :sigtau, :ell)
    all(haskey(vals, k) for k in needed) || return nothing
    return (; P=Int(vals[:P]), kappa=Float64(vals[:kappa]), tau0=Float64(vals[:tau0]),
             sigtau=Float64(vals[:sigtau]), ell=Float64(vals[:ell]))
end

function range_centered(c, halfwidth, step)
    lo = c - halfwidth
    hi = c + halfwidth
    # collect with rounding to avoid values like 0.048499999999
    return [round(x; digits=10) for x in collect(lo:step:hi)]
end

function robustness_center()
    if ROBUSTNESS_USE_BEST_FILE
        b = parse_best_so_far(ROBUSTNESS_BEST_FILE)
        if b !== nothing
            return b
        end
    end
    return (; P=ROBUST_P, kappa=ROBUST_KAPPA, tau0=ROBUST_TAU0, sigtau=ROBUST_SIGTAU, ell=ROBUST_ELL)
end

function robustness_grid_from_center(c)
    return (;
        P_values       = QUICK ? [c.P] : collect(max(1, c.P - ROBUST_DP):(c.P + ROBUST_DP)),
        kappa_values   = QUICK ? [c.kappa] : range_centered(c.kappa, ROBUST_DKAPPA, ROBUST_STEP_KAPPA),
        tau0_values    = QUICK ? [c.tau0] : range_centered(c.tau0, ROBUST_DTAU0, ROBUST_STEP_TAU0),
        sigtau_values  = QUICK ? [c.sigtau] : range_centered(c.sigtau, ROBUST_DSIGTAU, ROBUST_STEP_SIGTAU),
        ell_values     = [c.ell],
        delay_seeds    = ROBUST_DELAY_SEEDS,
        init_seeds     = ROBUST_INIT_SEEDS,
    )
end

function active_grid()
    if RUN_MODE == :robustness
        return robustness_grid_from_center(robustness_center())
    elseif RUN_MODE == :paper
        # Paper mode does not use this grid; return search grid only for harmless printing.
        return SEARCH_GRID
    else
        return SEARCH_GRID
    end
end
