# Run ordinary parameter search
# ENV["CHIMERA_MODE"] = "search"
include(joinpath(@__DIR__, "..", "src", "SarfatiPelegChimera.jl"))
SarfatiPelegChimera.main()
