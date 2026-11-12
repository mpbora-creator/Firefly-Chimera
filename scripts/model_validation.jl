# Standalone numerical tests of the Section-5 reduced delay-field model.
# Quick smoke test:
#   CHIMERA_QUICK=1 julia --project=. scripts/model_validation.jl
# Full validation:
#   julia --project=. scripts/model_validation.jl

ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl"))
using .SarfatiPelegChimera

result = SarfatiPelegChimera.run_model_validation()
fit = result.fit.selected
cfit = result.circular_fit.selected
println("Validation output: ",result.outdir)
println("Event-map test RMSE: ",fit.test.rmse)
println("Event-map test R²: ",fit.test.r2)
println("Selected held-out model: ",fit.name)
println("Features: ",SarfatiPelegChimera.REDUCED_FEATURES)
println("Coefficients: ",fit.coefficients)
println("Selected circular model: ",cfit.name)
println("Circular angular RMSE: ",cfit.metrics.angular_rmse)
println("Circular mean resultant: ",cfit.metrics.mean_resultant)
println("Transition metrics: ",result.circular_fit.probabilistic.transition)
