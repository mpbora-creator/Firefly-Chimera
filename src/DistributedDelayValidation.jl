# Validation of the distributed-delay stability analysis (Sections 9--14).
# Uses the actual link delays and weights of the numerical ring.

Base.@kwdef struct DistributedDelayConfig
    diffusion::Float64 = 0.0
    max_mode::Int = QUICK ? 8 : 50
    root_real_min::Float64 = -2.0
    root_real_max::Float64 = 0.5
    root_imag_max::Float64 = 8.0
    root_real_seeds::Int = QUICK ? 5 : 9
    root_imag_seeds::Int = QUICK ? 17 : 49
    hopf_omega_max::Float64 = 12.0
    hopf_samples::Int = QUICK ? 500 : 3000
    phase_mixing_samples::Int = QUICK ? 120 : 500
    residual_tolerance::Float64 = 1e-7
end

function distributed_transform(coupling,mode,lambda)
    N=length(coupling.neighbors); value=0.0+0im
    for i in 1:N,q in eachindex(coupling.neighbors[i])
        j=coupling.neighbors[i][q]; delta=j-i
        delta>N÷2&&(delta-=N); delta< -N÷2&&(delta+=N)
        value+=coupling.weights[i][q]*exp(-2pi*im*mode*delta/N-lambda*coupling.taus[i][q])
    end
    return value/N
end

function distributed_transform_derivative(coupling,mode,lambda)
    N=length(coupling.neighbors); value=0.0+0im
    for i in 1:N,q in eachindex(coupling.neighbors[i])
        j=coupling.neighbors[i][q]; delta=j-i
        delta>N÷2&&(delta-=N); delta< -N÷2&&(delta+=N); tau=coupling.taus[i][q]
        value-=coupling.weights[i][q]*tau*exp(-2pi*im*mode*delta/N-lambda*tau)
    end
    return value/N
end

function delay_moment(coupling,mode,order)
    N=length(coupling.neighbors); value=0.0+0im
    for i in 1:N,q in eachindex(coupling.neighbors[i])
        j=coupling.neighbors[i][q]; delta=j-i
        delta>N÷2&&(delta-=N); delta< -N÷2&&(delta+=N)
        value+=coupling.weights[i][q]*coupling.taus[i][q]^order*exp(-2pi*im*mode*delta/N)
    end
    return value/N
end

function delay_architecture_rows(coupling,modes)
    rows=NamedTuple[]
    for m in modes
        K=delay_moment(coupling,m,0); T1=delay_moment(coupling,m,1)
        T2=delay_moment(coupling,m,2); T3=delay_moment(coupling,m,3)
        push!(rows,(;mode=m,K_real=real(K),K_imag=imag(K),T1_real=real(T1),T1_imag=imag(T1),
                    T2_real=real(T2),T2_imag=imag(T2),T3_real=real(T3),T3_imag=imag(T3)))
    end
    return rows
end

function characteristic_residual(lambda,qstar,D,k,Lambda,coupling,mode)
    return lambda-qstar+D*k^2-Lambda*distributed_transform(coupling,mode,lambda)
end

function characteristic_newton(seed,qstar,D,k,Lambda,coupling,mode;tol=1e-9)
    lambda=ComplexF64(seed)
    for _ in 1:100
        F=characteristic_residual(lambda,qstar,D,k,Lambda,coupling,mode)
        deriv=1-Lambda*distributed_transform_derivative(coupling,mode,lambda)
        abs(deriv)<1e-12&&return nothing
        step=F/deriv; lambda-=step
        (!isfinite(real(lambda))||!isfinite(imag(lambda)))&&return nothing
        abs(step)<tol&&break
    end
    return lambda
end

function characteristic_roots(qstar,D,k,Lambda,coupling,mode,dcfg)
    roots=ComplexF64[]
    reals=range(dcfg.root_real_min,dcfg.root_real_max,length=dcfg.root_real_seeds)
    imags=range(-dcfg.root_imag_max,dcfg.root_imag_max,length=dcfg.root_imag_seeds)
    for x in reals,y in imags
        root=characteristic_newton(x+im*y,qstar,D,k,Lambda,coupling,mode)
        root===nothing&&continue
        abs(characteristic_residual(root,qstar,D,k,Lambda,coupling,mode))<=dcfg.residual_tolerance||continue
        all(abs(root-r)>1e-5 for r in roots)&&push!(roots,root)
    end
    sort!(roots,by=real,rev=true); return roots
end

function constant_delay_roots(qstar,D,k,Lambda,Khat,taumean,dcfg)
    roots=ComplexF64[]
    function newton(seed)
        z=ComplexF64(seed)
        for _ in 1:100
            F=z-qstar+D*k^2-Lambda*Khat*exp(-z*taumean)
            der=1+Lambda*Khat*taumean*exp(-z*taumean); abs(der)<1e-12&&return nothing
            step=F/der; z-=step; abs(step)<1e-9&&break
        end
        return z
    end
    for x in range(dcfg.root_real_min,dcfg.root_real_max,length=dcfg.root_real_seeds),y in range(-dcfg.root_imag_max,dcfg.root_imag_max,length=dcfg.root_imag_seeds)
        z=newton(x+im*y); z===nothing&&continue
        abs(z-qstar+D*k^2-Lambda*Khat*exp(-z*taumean))<dcfg.residual_tolerance||continue
        all(abs(z-r)>1e-5 for r in roots)&&push!(roots,z)
    end
    sort!(roots,by=real,rev=true); return roots
end

function moment_approximations(qstar,D,k,Lambda,coupling,mode)
    K=delay_moment(coupling,mode,0); T1=delay_moment(coupling,mode,1); T2=delay_moment(coupling,mode,2)
    numerator=qstar-D*k^2+Lambda*K
    first=numerator/(1+Lambda*T1)
    A=Lambda*T2/2; B=-(1+Lambda*T1); C=numerator
    if abs(A)<1e-14
        second=ComplexF64[first]
    else
        disc=B^2-4A*C; second=ComplexF64[(-B+sqrt(disc))/(2A),(-B-sqrt(disc))/(2A)]
    end
    return (;first,second,K,T1,T2)
end

function hopf_candidates(qstar,D,k,coupling,mode,dcfg)
    omegas=collect(range(1e-5,dcfg.hopf_omega_max,length=dcfg.hopf_samples)); residual=fill(NaN,length(omegas)); gains=fill(NaN,length(omegas))
    for i in eachindex(omegas)
        H=distributed_transform(coupling,mode,im*omegas[i]); abs(real(H))<1e-10&&continue
        L=-(qstar-D*k^2)/real(H)
        if L>0; gains[i]=L; residual[i]=omegas[i]-L*imag(H); end
    end
    rows=NamedTuple[]
    for i in 1:length(omegas)-1
        isfinite(residual[i])&&isfinite(residual[i+1])&&residual[i]*residual[i+1]<=0||continue
        lo,hi=omegas[i],omegas[i+1]
        for _ in 1:50
            mid=(lo+hi)/2; H=distributed_transform(coupling,mode,im*mid); L=-(qstar-D*k^2)/real(H); fm=mid-L*imag(H)
            Hlo=distributed_transform(coupling,mode,im*lo); Llo=-(qstar-D*k^2)/real(Hlo); flo=lo-Llo*imag(Hlo)
            flo*fm<=0 ? (hi=mid) : (lo=mid)
        end
        w=(lo+hi)/2; H=distributed_transform(coupling,mode,im*w); L=-(qstar-D*k^2)/real(H)
        L>0&&push!(rows,(;mode,omega=w,Lambda_critical=L,H_real=real(H),H_imag=imag(H)))
    end
    return rows
end

function phase_mixing_rows(coupling,modes,dcfg)
    rows=NamedTuple[]
    for m in modes,w in range(0,dcfg.hopf_omega_max,length=dcfg.phase_mixing_samples)
        H=distributed_transform(coupling,m,im*w); K=delay_moment(coupling,m,0)
        tau0=real(delay_moment(coupling,0,1)); Hconst=K*exp(-im*w*tau0)
        push!(rows,(;mode=m,omega=w,magnitude=abs(H),phase=angle(H),constant_magnitude=abs(Hconst),constant_phase=angle(Hconst)))
    end
    return rows
end

function distributed_stability_rows(equilibria,params,coupling,N,dcfg)
    rows=NamedTuple[]
    for eq in equilibria
        qstar=eq.delay==0 ? -params.a : -params.a+3*params.b*eq.delay^2-5*params.c*eq.delay^4
        for m in 0:min(dcfg.max_mode,N÷2)
            k=2pi*m/N; roots=characteristic_roots(qstar,dcfg.diffusion,k,params.Lambda,coupling,m,dcfg)
            K=delay_moment(coupling,m,0); tau0=real(delay_moment(coupling,0,1)); croots=constant_delay_roots(qstar,dcfg.diffusion,k,params.Lambda,K,tau0,dcfg)
            approx=moment_approximations(qstar,dcfg.diffusion,k,params.Lambda,coupling,m)
            slow_second=approx.second[argmin(abs.(approx.second.-approx.first))]
            for (branch,r) in enumerate(roots)
                push!(rows,(;equilibrium=eq.name,delay=eq.delay,mode=m,branch,
                            growth_rate=real(r),frequency=imag(r),dominant=branch==1,
                            constant_growth=isempty(croots) ? NaN : real(croots[1]),
                            first_moment_growth=real(approx.first),
                            second_moment_growth=real(slow_second)))
            end
        end
    end
    return rows
end

function plot_distributed_delay(architecture,phase,stability,hopf,outpath)
    p1=plot(;xlabel="frequency ω",ylabel="|Hₖ(iω)|",title="Distributed-delay phase mixing")
    for m in sort(unique(getproperty.(phase,:mode)))
        r=filter(x->x.mode==m,phase); plot!(p1,getproperty.(r,:omega),getproperty.(r,:magnitude);label="mode $m")
    end
    dominant=filter(x->x.dominant,stability); p2=plot(;xlabel="ring mode",ylabel="dominant growth rate",title="Full distributed-delay roots")
    for name in unique(getproperty.(dominant,:equilibrium))
        r=filter(x->x.equilibrium==name,dominant); plot!(p2,getproperty.(r,:mode),getproperty.(r,:growth_rate);marker=:circle,label=string(name))
    end
    chorus=filter(x->x.dominant&&x.equilibrium=="chorus",stability)
    p3=plot(getproperty.(chorus,:mode),getproperty.(chorus,:growth_rate);marker=:circle,label="full",xlabel="ring mode",ylabel="growth rate",title="Approximation comparison")
    plot!(p3,getproperty.(chorus,:mode),getproperty.(chorus,:constant_growth);label="constant delay",ls=:dash)
    plot!(p3,getproperty.(chorus,:mode),getproperty.(chorus,:first_moment_growth);label="first moment",ls=:dot)
    p4=scatter(getproperty.(hopf,:mode),getproperty.(hopf,:Lambda_critical);marker_z=getproperty.(hopf,:omega),colorbar_title="ω",xlabel="ring mode",ylabel="critical Λ",title="Distributed-delay Hopf candidates",label=false)
    plt=plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm,bottom_margin=12Plots.mm)
    savefig(plt,outpath); return outpath
end

function run_distributed_delay_validation(;cfg=CFG,case=analytical_validation_case(),
        outdir=joinpath(FIGURES_DIR,"distributed_delay_validation"),dcfg=DistributedDelayConfig())
    mkpath(outdir); lc=build_limit_cycle(cfg); baseline=simulate_case(cfg,case,lc;keep_full=true)
    obs=delay_relaxation_observables(baseline)
    # Calibrate the additive model used in the attached analysis from a full kappa scan.
    kappas=QUICK ? [case.kappa-.005,case.kappa,case.kappa+.005] : collect(case.kappa-.01:.0025:case.kappa+.01)
    scan=NamedTuple[]
    for kap in kappas
        ccase=merge(case,(;kappa=Float64(kap))); r=simulate_case(cfg,ccase,lc;keep_full=true); o=delay_relaxation_observables(r)
        G0=isfinite(o.lambda0) ? real(distributed_transform(r.coupling,0,o.lambda0)) : NaN
        push!(scan,(;kappa=Float64(kap),lambda0=o.lambda0,G0))
    end
    additive=constrained_susceptibility_fit(scan,:additive)
    params=infer_delayed_coefficients(obs,baseline.coupling,case.kappa,additive)
    equilibria=[(;name="chorus",delay=0.0)]
    if isfinite(obs.dchar); push!(equilibria,(;name="character_plus",delay=obs.dchar)); push!(equilibria,(;name="character_minus",delay=-obs.dchar)); end
    modes=collect(0:min(dcfg.max_mode,cfg.N÷2)); architecture=delay_architecture_rows(baseline.coupling,modes)
    stability=distributed_stability_rows(equilibria,params,baseline.coupling,cfg.N,dcfg)
    hopf=reduce(vcat,[hopf_candidates(-params.a,dcfg.diffusion,2pi*m/cfg.N,baseline.coupling,m,dcfg) for m in modes];init=NamedTuple[])
    phase=phase_mixing_rows(baseline.coupling,unique([0,min(case.P,maximum(modes))]),dcfg)
    write_namedtuples_csv(joinpath(outdir,"delay_architecture.csv"),architecture)
    write_namedtuples_csv(joinpath(outdir,"characteristic_roots.csv"),stability)
    write_namedtuples_csv(joinpath(outdir,"hopf_candidates.csv"),hopf)
    write_namedtuples_csv(joinpath(outdir,"phase_mixing.csv"),phase)
    chorusroots=filter(x->x.equilibrium=="chorus",stability)
    summary=[(;a=params.a,Lambda=params.Lambda,b=params.b,c=params.c,physical=params.physical,
               observed_chorus_rate=obs.lambda0,predicted_chorus_rate=isempty(chorusroots) ? NaN : maximum(getproperty.(chorusroots,:growth_rate)),
               additive_test_rmse=additive.test.rmse,additive_test_r2=additive.test.r2,diffusion=dcfg.diffusion)]
    write_namedtuples_csv(joinpath(outdir,"distributed_delay_summary.csv"),summary)
    !isempty(stability)&&!isempty(hopf)&&plot_distributed_delay(architecture,phase,stability,hopf,joinpath(outdir,"distributed_delay_validation.png"))
    return (;baseline,obs,additive,params,architecture,stability,hopf,phase,outdir)
end
