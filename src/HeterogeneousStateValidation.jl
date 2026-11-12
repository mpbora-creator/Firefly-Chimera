# Heterogeneous stationary states of the reduced distributed-delay graph model.

Base.@kwdef struct HeterogeneousConfig
    diffusion::Float64 = 0.0
    max_newton_iterations::Int = 80
    residual_tolerance::Float64 = 1e-9
    damping::Float64 = 0.5
    lm_iterations::Int = 120
    lm_initial_regularization::Float64 = 1e-4
end

function read_single_csv_row(path)
    lines=readlines(path); length(lines)>=2||error("missing data row in $path")
    names=Symbol.(split(lines[1],',')); values=split(lines[2],',')
    return Dict(names[i]=>values[i] for i in eachindex(names))
end

function distributed_parameter_file()
    joinpath(FIGURES_DIR,"distributed_delay_validation","distributed_delay_summary.csv")
end

function load_distributed_parameters(path=distributed_parameter_file())
    row=read_single_csv_row(path)
    getf(k)=parse(Float64,row[k])
    return (;a=getf(:a),Lambda=getf(:Lambda),b=getf(:b),c=getf(:c))
end

function graph_weight_matrix(coupling)
    N=length(coupling.neighbors); W=zeros(Float64,N,N)
    for i in 1:N,q in eachindex(coupling.neighbors[i]); W[i,coupling.neighbors[i][q]]+=coupling.weights[i][q]; end
    return W
end

function periodic_laplacian(N)
    L=zeros(Float64,N,N)
    for i in 1:N; L[i,i]=-2; L[i,mod1(i-1,N)]=1; L[i,mod1(i+1,N)]=1; end
    return L
end

function stationary_residual(d,p,W,L,D)
    return -p.a.*d .+ p.b.*d.^3 .- p.c.*d.^5 .+ D.*(L*d) .+ p.Lambda.*(W*d)
end

function stationary_jacobian(d,p,W,L,D)
    return Diagonal(-p.a .+3*p.b.*d.^2 .-5*p.c.*d.^4)+D.*L+p.Lambda.*W
end

function solve_heterogeneous_state(initial,p,coupling;hcfg=HeterogeneousConfig())
    d=Float64.(initial); W=graph_weight_matrix(coupling); L=periodic_laplacian(length(d))
    history=Float64[]
    for _ in 1:hcfg.max_newton_iterations
        F=stationary_residual(d,p,W,L,hcfg.diffusion); push!(history,norm(F,Inf))
        last(history)<hcfg.residual_tolerance&&break
        J=stationary_jacobian(d,p,W,L,hcfg.diffusion); step=try J\F catch; pinv(Matrix(J))*F end
        alpha=hcfg.damping; old=norm(F)
        while alpha>1e-5
            trial=d.-alpha.*step
            if norm(stationary_residual(trial,p,W,L,hcfg.diffusion))<old; d=trial; break; end
            alpha*=0.5
        end
        alpha<=1e-5&&break
    end
    # Levenberg-Marquardt refinement minimizes the residual even when Newton's
    # Jacobian is nearly singular at an interface.
    reg=hcfg.lm_initial_regularization
    for _ in 1:hcfg.lm_iterations
        F=stationary_residual(d,p,W,L,hcfg.diffusion); norm(F,Inf)<hcfg.residual_tolerance&&break
        J=Matrix(stationary_jacobian(d,p,W,L,hcfg.diffusion))
        step=(J'J+reg*I)\(J'F); trial=d.-step
        if norm(stationary_residual(trial,p,W,L,hcfg.diffusion))<norm(F)
            d=trial; reg=max(reg/3,1e-12); push!(history,norm(stationary_residual(d,p,W,L,hcfg.diffusion),Inf))
        else
            reg=min(reg*10,1e12)
        end
    end
    residual=norm(stationary_residual(d,p,W,L,hcfg.diffusion),Inf)
    return (;state=d,residual,converged=residual<hcfg.residual_tolerance,history,W,L)
end

function solve_heterogeneous_candidates(initial,p,coupling;hcfg=HeterogeneousConfig())
    mu=p.Lambda-p.a; disc=p.b^2+4p.c*mu
    dchar=disc>=0 ? sqrt(max(0,(p.b+sqrt(disc))/(2p.c))) : maximum(abs.(initial))
    quantized=[abs(x)<CHARACTER_TRACK_MIN_DELAY ? 0.0 : sign(x)*dchar for x in initial]
    starts=[initial,0.25.*initial,0.5.*initial,quantized]
    solutions=[solve_heterogeneous_state(s,p,coupling;hcfg) for s in starts]
    return solutions[argmin(getproperty.(solutions,:residual))],solutions
end

function numerical_delay_profile(result)
    D=result.delay_tracks.delay_matrix; profile=zeros(size(D,1))
    for i in axes(D,1)
        vals=Float64[x for x in view(D,i,:) if isfinite(x)]; profile[i]=isempty(vals) ? 0.0 : mean(vals)
    end
    return profile
end

function heterogeneous_profile_metrics(numerical,reduced,threshold)
    corr=std(numerical)>0&&std(reduced)>0 ? cor(numerical,reduced) : NaN
    rmse=sqrt(mean((numerical.-reduced).^2))
    return (;correlation=corr,rmse,numerical_character_fraction=mean(abs.(numerical).>threshold),
             reduced_character_fraction=mean(abs.(reduced).>threshold))
end

function run_heterogeneous_state_validation(;cfg=CFG,case=analytical_validation_case(),
        outdir=joinpath(FIGURES_DIR,"heterogeneous_state_validation"),hcfg=HeterogeneousConfig())
    mkpath(outdir); p=load_distributed_parameters(); lc=build_limit_cycle(cfg)
    result=simulate_case(cfg,case,lc;keep_full=true); numerical=numerical_delay_profile(result)
    solved,candidates=solve_heterogeneous_candidates(numerical,p,result.coupling;hcfg)
    threshold=max(CFG.chorus_window,CHARACTER_TRACK_MIN_DELAY)
    metrics=heterogeneous_profile_metrics(numerical,solved.state,threshold)
    rows=[(;oscillator=i,numerical_delay=numerical[i],reduced_delay=solved.state[i]) for i in eachindex(numerical)]
    write_namedtuples_csv(joinpath(outdir,"heterogeneous_profile.csv"),rows)
    write_namedtuples_csv(joinpath(outdir,"heterogeneous_metrics.csv"),[(;metrics...,residual=solved.residual,converged=solved.converged,diffusion=hcfg.diffusion)])
    p1=plot(eachindex(numerical),numerical;label="numerical",xlabel="oscillator",ylabel="mean delay",title="Heterogeneous state")
    plot!(p1,eachindex(solved.state),solved.state;label="reduced stationary",lw=2)
    p2=scatter(numerical,solved.state;label=false,xlabel="numerical delay",ylabel="reduced delay",title="Profile agreement")
    lim=collect(extrema(vcat(numerical,solved.state))); plot!(p2,lim,lim;ls=:dash,label="identity")
    p3=plot(solved.history;yscale=:log10,label=false,xlabel="Newton iteration",ylabel="residual",title="Stationary convergence")
    p4=histogram(numerical;bins=25,alpha=.5,label="numerical",xlabel="delay",ylabel="count",title="Delay-state distribution")
    histogram!(p4,solved.state;bins=25,alpha=.5,label="reduced")
    savefig(plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm),joinpath(outdir,"heterogeneous_state_validation.png"))
    return (;result,p,solved,candidates,metrics,outdir)
end
