# Sarfati–Peleg Chimera Framework
# Shared loader for search, robustness, and paper-figure workflows.

module SarfatiPelegChimera

include("Config.jl")
include("Utilities.jl")
include("LimitCycle.jl")
include("Coupling.jl")
include("Simulation.jl")
include("Diagnostics.jl")
include("ModelValidation.jl")
include("AnalyticalValidation.jl")
include("RevisedAnalyticalValidation.jl")
include("DistributedDelayValidation.jl")
include("HeterogeneousStateValidation.jl")
include("ArgumentPrincipleValidation.jl")
include("ReductionRobustnessValidation.jl")
include("TimeDependentReducedValidation.jl")
include("LoggingResume.jl")
include("Plotting.jl")
include("Robustness.jl")
include("PaperFigures.jl")
include("Main.jl")

end # module
