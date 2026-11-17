#!/usr/bin/env julia

using CSV
using DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUT = joinpath(ROOT, "results", "natural_chimera_results.csv")
const OUTDIR = joinpath(ROOT, "figures", "role_switching_robustness")

mkpath(OUTDIR)
d = CSV.read(INPUT, DataFrame)

# Keep every non-scanned quantity at the representative Fig. 10 values.
slice = d[(d.sigtau .== 0.062) .&
          (d.delay_seed .== 29) .&
          (d.init_seed .== 202), :]

sort!(slice, [:tau0, :kappa])
out = select(slice, :kappa, :tau0, :mean_character_probability,
             :switching_oscillator_count, :persistent_chorus_count,
             :persistent_character_count, :mean_role_switch_count)
out.switching_fraction = out.switching_oscillator_count ./ 100
CSV.write(joinpath(OUTDIR, "role_switching_phase_map.dat"), out; delim=' ')

println("Saved ", nrow(out), " parameter points to ", OUTDIR)
