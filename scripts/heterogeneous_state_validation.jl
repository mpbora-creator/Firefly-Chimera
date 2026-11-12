ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl")); using .SarfatiPelegChimera
r=SarfatiPelegChimera.run_heterogeneous_state_validation(); println("Output: ",r.outdir," metrics=",r.metrics," residual=",r.solved.residual)
