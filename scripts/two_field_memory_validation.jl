# Standalone two-field delay-memory validation. Existing source is loaded read-only.
ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl"))
include(joinpath(@__DIR__,"..","src","TwoFieldMemoryModelValidation.jl"))
using .TwoFieldMemoryModelValidation
r=TwoFieldMemoryModelValidation.run()
println("Two-field validation output: ",r.outdir)
println("Selected parameters: ",r.fit.coefficients," beta=",r.fit.beta," gamma=",r.fit.gamma," tau_m=",r.fit.tau_m)
println("Held-out metrics: ",r.fit.test_metrics)
println("Equilibria: ",r.equilibria)
println("Spectrum roots: ",length(r.spectrum))
