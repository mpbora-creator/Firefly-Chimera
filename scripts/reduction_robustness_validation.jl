ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl")); using .SarfatiPelegChimera
r=SarfatiPelegChimera.run_reduction_robustness_validation(); println("Output: ",r.outdir," cases=",length(r.rows))
