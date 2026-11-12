# Standalone settle-kick-observe basin-boundary validation.
ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl"))
include(joinpath(@__DIR__,"..","src","TwoFieldBasinBoundaryValidation.jl"))
using .TwoFieldBasinBoundaryValidation
r=TwoFieldBasinBoundaryValidation.run()
println("Basin validation output: ",r.outdir)
println("Analytical unstable boundary: ",r.boundary)
println("Numerical thresholds: ",r.thresholds)
