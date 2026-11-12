# Regenerate full diagnostics for top completed cases
# ENV["CHIMERA_MODE"] = "paper"
include(joinpath(@__DIR__, "..", "src", "SarfatiPelegChimera.jl"))
SarfatiPelegChimera.main()
