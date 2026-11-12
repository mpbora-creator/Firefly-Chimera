# Revised delayed/nonlocal/periodic analytical validation.
# CHIMERA_QUICK=1 julia --project=. scripts/revised_analytical_validation.jl
# julia --project=. scripts/revised_analytical_validation.jl

ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl"))
using .SarfatiPelegChimera

r=SarfatiPelegChimera.run_revised_analytical_validation()
println("Revised analytical output: ",r.outdir)
println("Selected coupling form: ",r.susceptibility.form)
println("Held-out susceptibility metrics: ",r.susceptibility.test)
println("Coupling comparison: ",[(f.form,f.test) for f in r.comparison.fits])
println("Delayed coefficients: ",r.coeff)
println("Periodic coefficients: ",r.periodic)
println("Periodic equilibria: ",r.eqs)
