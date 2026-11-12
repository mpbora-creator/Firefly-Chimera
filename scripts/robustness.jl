# Run local robustness/sensitivity study around best-so-far or manual centre
# ENV["CHIMERA_MODE"] = "robustness"
include(joinpath(@__DIR__, "..", "src", "SarfatiPelegChimera.jl"))
SarfatiPelegChimera.main()
