# Full distributed-delay stability validation (Sections 9--14).
# CHIMERA_QUICK=1 julia --project=. scripts/distributed_delay_validation.jl
# julia --project=. scripts/distributed_delay_validation.jl

ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl"))
using .SarfatiPelegChimera

r=SarfatiPelegChimera.run_distributed_delay_validation()
println("Distributed-delay output: ",r.outdir)
println("Additive held-out metrics: ",r.additive.test)
println("Parameters: ",r.params)
println("Roots found: ",length(r.stability),"  Hopf candidates: ",length(r.hopf))
