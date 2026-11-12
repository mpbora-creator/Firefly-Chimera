# Multi-seed robustness of numerical-to-reduced heterogeneous profile agreement.

Base.@kwdef struct ReductionRobustnessConfig
    delay_seeds::Vector{Int} = QUICK ? [11,29] : [11,17,29,37,43]
    init_seeds::Vector{Int} = QUICK ? [101] : [101,202,303]
end

function run_reduction_robustness_validation(;cfg=CFG,basecase=analytical_validation_case(),
        outdir=joinpath(FIGURES_DIR,"reduction_robustness_validation"),rcfg=ReductionRobustnessConfig(),hcfg=HeterogeneousConfig())
    mkpath(outdir); p=load_distributed_parameters(); lc=build_limit_cycle(cfg); rows=NamedTuple[]
    for ds in rcfg.delay_seeds,is in rcfg.init_seeds
        case=merge(basecase,(;delay_seed=ds,init_seed=is)); result=simulate_case(cfg,case,lc;keep_full=true)
        numerical=numerical_delay_profile(result); solved=solve_heterogeneous_state(numerical,p,result.coupling;hcfg)
        met=heterogeneous_profile_metrics(numerical,solved.state,max(cfg.chorus_window,CHARACTER_TRACK_MIN_DELAY))
        push!(rows,(;delay_seed=ds,init_seed=is,met...,residual=solved.residual,converged=solved.converged,
                    score=result.summary.score,stable_character_fraction=result.summary.stable_character_fraction))
    end
    write_namedtuples_csv(joinpath(outdir,"reduction_robustness.csv"),rows)
    p1=scatter(getproperty.(rows,:score),getproperty.(rows,:correlation);marker_z=getproperty.(rows,:stable_character_fraction),
        xlabel="numerical score",ylabel="profile correlation",colorbar_title="stable character fraction",label=false,title="Reduction robustness")
    p2=bar(eachindex(rows),getproperty.(rows,:rmse);legend=false,xlabel="seed pair",ylabel="profile RMSE",title="Profile error")
    savefig(plot(p1,p2;layout=(1,2),size=(1000,420),margin=7Plots.mm),joinpath(outdir,"reduction_robustness_validation.png"))
    return (;rows,outdir)
end
