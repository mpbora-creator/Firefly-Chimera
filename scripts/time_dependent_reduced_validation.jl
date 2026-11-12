ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl")); using .SarfatiPelegChimera
r=SarfatiPelegChimera.run_time_dependent_reduced_validation(); println("Output: ",r.outdir," metrics=",r.metrics," period=",r.period)
