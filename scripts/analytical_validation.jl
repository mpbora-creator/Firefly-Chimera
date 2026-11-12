# Forward validation of the Section-5 analytical tri-stability model.
# CHIMERA_QUICK=1 julia --project=. scripts/analytical_validation.jl
# julia --project=. scripts/analytical_validation.jl

ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl"))
using .SarfatiPelegChimera

r=SarfatiPelegChimera.run_analytical_validation()
println("Analytical validation output: ",r.outdir)
println("Measured observables: ",r.obs)
println("Inferred coefficients: ",r.coeff)
println("Bifurcation predictions: ",r.bif)
println("Measured basin threshold: ",r.basin.measured_threshold)
println("Correlation length: ",r.covfit.xi,"  effective D: ",r.D_eff)
