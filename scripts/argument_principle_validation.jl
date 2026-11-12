ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl")); using .SarfatiPelegChimera
r=SarfatiPelegChimera.run_argument_principle_validation(); println("Output: ",r.outdir," certificate=",r.certificate===nothing ? "refused" : r.certificate.unstable_roots)
