# Time-dependent reduced distributed-delay field. This tests whether the
# numerical chimera corresponds to a periodic/statistically stationary profile
# rather than a stationary solution.

Base.@kwdef struct TimeDependentReducedConfig
    dt::Float64 = 0.02
    t_end::Float64 = QUICK ? 80.0 : 400.0
    transient_fraction::Float64 = 0.5
    sample_dt::Float64 = 0.2
    diffusion::Float64 = 0.0
    perturbation::Float64 = 1e-7
end

function reduced_delayed_value(history,current,delay,dt,j)
    pos=current-delay/dt
    pos<=1&&return history[j,1]
    k0=clamp(floor(Int,pos),1,current); k1=clamp(k0+1,1,current); alpha=pos-k0
    return (1-alpha)*history[j,k0]+alpha*history[j,k1]
end

function integrate_reduced_delay(initial,p,coupling;tcfg=TimeDependentReducedConfig())
    N=length(initial); steps=round(Int,tcfg.t_end/tcfg.dt)+1; every=max(1,round(Int,tcfg.sample_dt/tcfg.dt))
    history=zeros(Float64,N,steps); history[:,1].=initial; d=copy(Float64.(initial)); L=periodic_laplacian(N)
    samples=Vector{Vector{Float64}}(); times=Float64[]
    for step in 1:steps-1
        rhs=-p.a.*d .+p.b.*d.^3 .-p.c.*d.^5 .+tcfg.diffusion.*(L*d)
        for i in 1:N
            delayed=sum(coupling.weights[i][q]*reduced_delayed_value(history,step,coupling.taus[i][q],tcfg.dt,coupling.neighbors[i][q]) for q in eachindex(coupling.neighbors[i]))
            rhs[i]+=p.Lambda*delayed
        end
        d .+=tcfg.dt.*rhs; history[:,step+1].=d
        if step%every==0; push!(samples,copy(d)); push!(times,step*tcfg.dt); end
    end
    return (;times,samples=isempty(samples) ? zeros(N,0) : reduce(hcat,samples),history)
end

function dominant_period(times,signal,start)
    ids=findall(times.>=start); length(ids)>=8||return NaN
    x=signal[ids].-mean(signal[ids]); maxlag=fld(length(x),2); maxlag>=2||return NaN
    ac=[dot(view(x,1:length(x)-l),view(x,1+l:length(x))) for l in 1:maxlag]
    peaks=[l for l in 2:maxlag-1 if ac[l]>=ac[l-1]&&ac[l]>=ac[l+1]]
    isempty(peaks) ? NaN : times[ids[min(first(peaks)+1,length(ids))]]-times[ids[1]]
end

function run_time_dependent_reduced_validation(;cfg=CFG,case=analytical_validation_case(),
        outdir=joinpath(FIGURES_DIR,"time_dependent_reduced_validation"),tcfg=TimeDependentReducedConfig())
    mkpath(outdir); p=load_distributed_parameters(); lc=build_limit_cycle(cfg); numerical=simulate_case(cfg,case,lc;keep_full=true)
    initial=numerical_delay_profile(numerical)
    mu=p.Lambda-p.a; disc=p.b^2+4p.c*mu
    dchar=disc>=0 ? sqrt(max(0,(p.b+sqrt(disc))/(2p.c))) : maximum(abs.(initial))
    # Event delays can differ by whole chorus periods. Start the reduced local
    # field on the physically equivalent principal branch of its normal form.
    reduced_initial=clamp.(initial,-2dchar,2dchar)
    reduced=integrate_reduced_delay(reduced_initial,p,numerical.coupling;tcfg)
    start=tcfg.transient_fraction*tcfg.t_end; ids=findall(reduced.times.>=start)
    meanprofile=vec(mean(reduced.samples[:,ids],dims=2)); stdprofile=vec(std(reduced.samples[:,ids],dims=2))
    metrics=heterogeneous_profile_metrics(initial,meanprofile,max(cfg.chorus_window,CHARACTER_TRACK_MIN_DELAY))
    globalmean=vec(mean(reduced.samples,dims=1)); period=dominant_period(reduced.times,globalmean,start)
    rows=[(;oscillator=i,numerical_delay=initial[i],reduced_mean=meanprofile[i],reduced_std=stdprofile[i]) for i in eachindex(initial)]
    write_namedtuples_csv(joinpath(outdir,"time_dependent_profile.csv"),rows)
    write_namedtuples_csv(joinpath(outdir,"time_dependent_metrics.csv"),[(;metrics...,dominant_period=period,
        mean_temporal_std=mean(stdprofile),max_temporal_std=maximum(stdprofile),diffusion=tcfg.diffusion)])
    p1=plot(eachindex(initial),initial;label="numerical",xlabel="oscillator",ylabel="delay",title="Time-averaged reduced profile")
    plot!(p1,eachindex(meanprofile),meanprofile;ribbon=stdprofile,label="reduced mean ± std")
    p2=heatmap(reduced.times,1:cfg.N,reduced.samples;xlabel="time",ylabel="oscillator",title="Reduced delay dynamics")
    p3=plot(reduced.times,globalmean;label=false,xlabel="time",ylabel="spatial mean delay",title="Temporal dynamics")
    p4=scatter(initial,meanprofile;label=false,xlabel="numerical delay",ylabel="reduced time mean",title="Profile agreement")
    savefig(plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm),joinpath(outdir,"time_dependent_reduced_validation.png"))
    return (;numerical,reduced,meanprofile,stdprofile,metrics,period,outdir)
end
