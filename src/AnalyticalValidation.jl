# Forward validation of the Section-5 tri-stable analytical model.
# Numerical observables calibrate a small set of coefficients; the analytical
# equations then predict basin thresholds, bifurcation points, and length scales.

Base.@kwdef struct AnalyticalValidationConfig
    near_chorus_window::Float64 = 0.35
    minimum_rate_samples::Int = 12
    covariance_max_distance::Int = 12
    basin_character_fraction::Float64 = 0.02
end

function analytical_validation_case(path=joinpath(RESULTS_DIR,"natural_chimera_best_so_far.txt"))
    isfile(path)||return SAMPLE_CASE
    vals=Dict{Symbol,Float64}()
    for line in readlines(path)
        occursin("=",line)||continue
        key,value=strip.(split(line,"=";limit=2))
        try vals[Symbol(key)]=parse(Float64,value) catch end
    end
    needed=(:P,:kappa,:tau0,:sigtau,:ell,:delay_seed,:init_seed)
    all(haskey(vals,k) for k in needed)||return SAMPLE_CASE
    return (;P=round(Int,vals[:P]),kappa=vals[:kappa],tau0=vals[:tau0],sigtau=vals[:sigtau],ell=vals[:ell],
             delay_seed=round(Int,vals[:delay_seed]),init_seed=round(Int,vals[:init_seed]))
end

function slope_through_origin(x,y)
    den=sum(abs2,x)
    return den>eps() ? dot(x,y)/den : NaN
end

function delay_relaxation_observables(result; acfg=AnalyticalValidationConfig())
    D=result.delay_tracks.delay_matrix; chorus=result.chorus
    N,C=size(D); chorus_x=Float64[]; chorus_y=Float64[]
    for n in 1:C-1, i in 1:N
        d0,d1=Float64(D[i,n]),Float64(D[i,n+1]); dt=chorus[n+1]-chorus[n]
        if isfinite(d0)&&isfinite(d1)&&dt>0&&abs(d0)<=acfg.near_chorus_window
            push!(chorus_x,d0); push!(chorus_y,(d1-d0)/dt)
        end
    end
    lambda0=length(chorus_x)>=acfg.minimum_rate_samples ? slope_through_origin(chorus_x,chorus_y) : NaN

    char_x=Float64[]; char_y=Float64[]; amplitudes=Float64[]
    ids=result.delay_tracks.stable_oscillators
    starts=result.delay_tracks.stable_start_cycles; stops=result.delay_tracks.stable_stop_cycles
    for q in eachindex(ids)
        i=ids[q]; lo=max(1,starts[q]); hi=min(C,stops[q]); hi-lo>=2 || continue
        vals=Float64[D[i,n] for n in lo:hi if isfinite(D[i,n])]
        isempty(vals)&&continue
        dstar=mean(vals); push!(amplitudes,abs(dstar))
        for n in lo:hi-1
            d0,d1=Float64(D[i,n]),Float64(D[i,n+1]); dt=chorus[n+1]-chorus[n]
            if isfinite(d0)&&isfinite(d1)&&dt>0
                delta=d0-dstar; push!(char_x,delta); push!(char_y,(d1-d0)/dt)
            end
        end
    end
    lambdachar=length(char_x)>=acfg.minimum_rate_samples ? slope_through_origin(char_x,char_y) : NaN
    dchar=isempty(amplitudes) ? NaN : median(amplitudes)
    return (;lambda0,lambdachar,dchar,chorus_samples=length(chorus_x),
             character_samples=length(char_x),character_count=length(amplitudes))
end

function infer_tristable_coefficients(obs)
    mu=obs.lambda0; z=obs.dchar^2; lp=obs.lambdachar
    if !(isfinite(mu)&&isfinite(z)&&z>0&&isfinite(lp))
        return (;mu,b=NaN,c=NaN,zplus=z,zminus=NaN,d_unstable=NaN,
                 mu_sn=NaN,d_sn=NaN,physically_consistent=false)
    end
    c=-(mu+lp/2)/z^2
    b=c*z-mu/z
    disc=b^2+4c*mu
    zminus=c>0&&disc>=0 ? (b-sqrt(disc))/(2c) : NaN
    dunstable=isfinite(zminus)&&zminus>0 ? sqrt(zminus) : NaN
    physical=b>0&&c>0&&mu<0&&disc>=0
    return (;mu,b,c,zplus=z,zminus,d_unstable=dunstable,
             mu_sn=c>0 ? -b^2/(4c) : NaN,d_sn=b>0&&c>0 ? sqrt(b/(2c)) : NaN,
             physically_consistent=physical)
end

function covariance_length_scale(result; acfg=AnalyticalValidationConfig())
    r=result.delay_cov_distance.distances; cov=result.delay_cov_distance.mean_cov
    maxr=min(acfg.covariance_max_distance,maximum(r)); ids=[i for i in eachindex(r) if 1<=r[i]<=maxr&&isfinite(cov[i])&&cov[i]>0]
    length(ids)>=3 || return (;xi=NaN,slope=NaN,distances=Float64[],values=Float64[])
    x=Float64.(r[ids]); y=log.(cov[ids]); A=hcat(ones(length(x)),x); beta=A\y
    xi=beta[2]<0 ? -1/beta[2] : NaN
    return (;xi,slope=beta[2],distances=x,values=cov[ids])
end

function coupling_relaxation_scan(cfg,basecase,lc,kappas;acfg=AnalyticalValidationConfig())
    rows=NamedTuple[]
    for kappa in kappas
        case=merge(basecase,(;kappa=Float64(kappa)))
        result=simulate_case(cfg,case,lc;keep_full=true)
        obs=delay_relaxation_observables(result;acfg)
        push!(rows,(;kappa=Float64(kappa),lambda0=obs.lambda0,dchar=obs.dchar,
                    lambdachar=obs.lambdachar,stable_characters=obs.character_count))
    end
    valid=filter(x->isfinite(x.lambda0),rows)
    if length(valid)>=2
        X=hcat(ones(length(valid)),getproperty.(valid,:kappa)); y=getproperty.(valid,:lambda0); beta=X\y
        intercept,gamma=beta; a=-intercept
    else
        intercept=gamma=a=NaN
    end
    return (;rows,a,gamma,intercept)
end

function analytical_bifurcation_predictions(coeff,calibration)
    a,gamma=calibration.a,calibration.gamma
    if !(isfinite(a)&&isfinite(gamma)&&abs(gamma)>eps()&&isfinite(coeff.mu_sn))
        return (;kappa_sn=NaN,kappa_pf=NaN,hysteresis_width=NaN)
    end
    kpf=a/gamma; ksn=(a+coeff.mu_sn)/gamma
    return (;kappa_sn=ksn,kappa_pf=kpf,hysteresis_width=kpf-ksn)
end

function basin_threshold_scan(cfg,basecase,lc,offsets;acfg=AnalyticalValidationConfig())
    rows=NamedTuple[]
    for offset in offsets
        phases=controlled_phase_profile(cfg.N,:positive,0.15,Float64(offset))
        result=simulate_case(cfg,basecase,lc;keep_full=true,initial_phases=phases)
        D=result.delay_tracks.delay_matrix; tail=D[:,max(1,size(D,2)-max(3,round(Int,.25size(D,2)))+1):end]
        vals=tail[isfinite.(tail)]; outcome=result.summary.stable_character_fraction>=acfg.basin_character_fraction
        T=length(result.chorus)>=2 ? median(diff(result.chorus)) : NaN
        push!(rows,(;phase_offset=Float64(offset),initial_delay=isfinite(T) ? abs(offset)*T/(2pi) : NaN,
                    character_outcome=outcome,stable_character_fraction=result.summary.stable_character_fraction,
                    tail_mean_abs_delay=isempty(vals) ? NaN : mean(abs.(vals))))
    end
    hit=findfirst(x->x.character_outcome,rows)
    measured=hit===nothing ? NaN : rows[hit].initial_delay
    return (;rows,measured_threshold=measured)
end

function plot_analytical_validation(obs,coeff,calibration,bif,basin,covfit,outpath)
    kr=getproperty.(calibration.rows,:kappa); lr=getproperty.(calibration.rows,:lambda0)
    p1=scatter(kr,lr;label="measured",xlabel="κ",ylabel="chorus relaxation rate",title="μ(κ) calibration")
    if isfinite(calibration.gamma)
        kk=range(minimum(kr),maximum(kr),length=100); plot!(p1,kk,calibration.intercept .+ calibration.gamma.*kk;label="linear fit",lw=2)
        isfinite(bif.kappa_sn)&&vline!(p1,[bif.kappa_sn,bif.kappa_pf];ls=:dash,label=["κSN" "κPF"])
    end
    delays=getproperty.(basin.rows,:initial_delay); outcomes=Float64.(getproperty.(basin.rows,:character_outcome))
    p2=scatter(delays,outcomes;ms=5,label="numerical",xlabel="initial delay kick",ylabel="character outcome",yticks=([0,1],["chorus","character"]),title="Basin boundary")
    isfinite(coeff.d_unstable)&&vline!(p2,[coeff.d_unstable];ls=:dash,lw=2,label="analytical dᵤ")
    ampnames=["stable character","unstable boundary","fold"]; ampvals=[obs.dchar,coeff.d_unstable,coeff.d_sn]
    p3=bar(1:3,ampvals;xticks=(1:3,ampnames),xrotation=20,legend=false,ylabel="delay amplitude",title="Analytical amplitudes")
    p4=plot(;xlabel="ring distance",ylabel="delay covariance",yscale=:log10,title="Correlation length")
    if !isempty(covfit.distances)
        scatter!(p4,covfit.distances,covfit.values;label="numerical")
        fitvals=exp.(log(covfit.values[1]).+(covfit.distances.-covfit.distances[1]).*covfit.slope)
        plot!(p4,covfit.distances,fitvals;label="ξ=$(round(covfit.xi,digits=2))",lw=2)
    end
    plt=plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm,bottom_margin=12Plots.mm)
    savefig(plt,outpath); return outpath
end

function run_analytical_validation(;cfg=CFG,case=analytical_validation_case(),
        kappas=QUICK ? [case.kappa-.005,case.kappa,case.kappa+.005] : collect(case.kappa-.01:.0025:case.kappa+.01),
        offsets=QUICK ? collect(0.0:.35:1.4) : collect(0.0:.15:1.8),
        outdir=joinpath(FIGURES_DIR,"analytical_validation"),acfg=AnalyticalValidationConfig())
    mkpath(outdir); lc=build_limit_cycle(cfg)
    baseline=simulate_case(cfg,case,lc;keep_full=true)
    obs=delay_relaxation_observables(baseline;acfg); coeff=infer_tristable_coefficients(obs)
    calibration=coupling_relaxation_scan(cfg,case,lc,kappas;acfg)
    bif=analytical_bifurcation_predictions(coeff,calibration)
    basin=basin_threshold_scan(cfg,case,lc,offsets;acfg)
    covfit=covariance_length_scale(baseline;acfg)
    D_eff=isfinite(covfit.xi)&&isfinite(obs.lambda0) ? covfit.xi^2*abs(obs.lambda0) : NaN
    summary=[(;mu=coeff.mu,b=coeff.b,c=coeff.c,d_character=obs.dchar,
              lambda_character=obs.lambdachar,d_unstable_predicted=coeff.d_unstable,
              d_unstable_measured=basin.measured_threshold,d_fold=coeff.d_sn,
              kappa_sn=bif.kappa_sn,kappa_pf=bif.kappa_pf,hysteresis_width=bif.hysteresis_width,
              correlation_length=covfit.xi,D_effective=D_eff,physically_consistent=coeff.physically_consistent)]
    write_namedtuples_csv(joinpath(outdir,"analytical_predictions.csv"),summary)
    write_namedtuples_csv(joinpath(outdir,"coupling_relaxation_scan.csv"),calibration.rows)
    write_namedtuples_csv(joinpath(outdir,"basin_threshold_scan.csv"),basin.rows)
    plot_analytical_validation(obs,coeff,calibration,bif,basin,covfit,joinpath(outdir,"analytical_validation.png"))
    return (;baseline,obs,coeff,calibration,bif,basin,covfit,D_eff,outdir)
end
