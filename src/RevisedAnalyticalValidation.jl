# Revised forward model: delayed eigenvalues, nonlinear susceptibility,
# periodic tri-stability, exact graph modes, and localized nucleation.

Base.@kwdef struct RevisedAnalyticalConfig
    susceptibility_degree::Int = 3
    nucleus_steps::Int = QUICK ? 1200 : 5000
    nucleus_dt::Float64 = 0.02
    nucleus_widths::Vector{Float64} = QUICK ? [0.05,0.15,0.30] : collect(0.03:0.04:0.35)
    nucleus_amplitude_factors::Vector{Float64} = QUICK ? [0.7,1.0,1.3] : collect(0.5:0.1:1.5)
end

function delay_transform(coupling,lambda)
    vals=ComplexF64[]
    for i in eachindex(coupling.neighbors),q in eachindex(coupling.neighbors[i])
        push!(vals,coupling.weights[i][q]*exp(-lambda*coupling.taus[i][q]))
    end
    # Each row has unit total weight, hence divide by the number of rows.
    return sum(vals)/length(coupling.neighbors)
end

function susceptibility_fit(scan_results,degree)
    valid=filter(x->isfinite(x.lambda0),scan_results)
    kappas=getproperty.(valid,:kappa); lambdas=getproperty.(valid,:lambda0)
    G=getproperty.(valid,:G0)
    X=hcat(ones(length(valid)),[G[i]*kappas[i]^p for i in eachindex(valid),p in 1:degree])
    beta=X\lambdas; predicted=X*beta
    den=sum((lambdas.-mean(lambdas)).^2); r2=den>0 ? 1-sum((lambdas.-predicted).^2)/den : NaN
    a=-beta[1]; gamma=beta[2:end]
    return (;a,gamma,r2,kappas,lambdas,predicted,G)
end

effective_gain(kappa,gamma)=sum(gamma[p]*kappa^p for p in eachindex(gamma))

function nelder_mead(objective,x0;step=0.25,maxiter=1200,tol=1e-9)
    n=length(x0); simplex=[copy(x0) for _ in 1:n+1]
    for j in 1:n; simplex[j+1][j]+=step; end
    values=objective.(simplex)
    for _ in 1:maxiter
        order=sortperm(values); simplex=simplex[order]; values=values[order]
        maximum(abs.(values.-values[1]))<tol&&break
        centroid=sum(simplex[1:n])/n
        reflected=centroid+(centroid-simplex[end]); fr=objective(reflected)
        if values[1]<=fr<values[n]
            simplex[end]=reflected; values[end]=fr
        elseif fr<values[1]
            expanded=centroid+2*(reflected-centroid); fe=objective(expanded)
            simplex[end],values[end]=fe<fr ? (expanded,fe) : (reflected,fr)
        else
            contracted=centroid+0.5*(simplex[end]-centroid); fc=objective(contracted)
            if fc<values[end]
                simplex[end]=contracted; values[end]=fc
            else
                for j in 2:n+1
                    simplex[j]=simplex[1]+0.5*(simplex[j]-simplex[1]); values[j]=objective(simplex[j])
                end
            end
        end
    end
    i=argmin(values); return simplex[i],values[i]
end

function constrained_gain(kappa,params,kmean,kscale)
    x=(kappa-kmean)/kscale
    return kappa*exp(clamp(params[2]+params[3]*x+params[4]*x^2,-30.0,30.0))
end

function constrained_susceptibility_fit(scan,form::Symbol;train_fraction=0.70,ridge=1e-3)
    form in (:additive,:difference)||error("coupling form must be additive or difference")
    valid=filter(x->isfinite(x.lambda0)&&isfinite(x.G0),scan)
    m=length(valid); m>=2||error("need at least two valid relaxation measurements")
    order=sortperm(getproperty.(valid,:kappa)); cut=m<5 ? m-1 : clamp(round(Int,train_fraction*m),3,m-2)
    train=order[1:cut]; test=order[cut+1:end]
    kappas=getproperty.(valid,:kappa); lambdas=getproperty.(valid,:lambda0); G=getproperty.(valid,:G0)
    kmean=mean(kappas[train]); rawscale=length(train)>=2 ? std(kappas[train]) : std(kappas)
    kscale=isfinite(rawscale)&&rawscale>eps() ? rawscale : 1.0
    function residual(i,p)
        a=exp(clamp(p[1],-30.0,30.0)); L=constrained_gain(kappas[i],p,kmean,kscale)
        return form==:additive ? lambdas[i]+a-L*G[i] : lambdas[i]+a+L-L*G[i]
    end
    objective(p)=mean(residual(i,p)^2 for i in train)+ridge*(p[3]^2+p[4]^2)
    p,_=nelder_mead(objective,[log(0.05),0.0,0.0,0.0];step=.4)
    a=exp(p[1]); prediction=[form==:additive ? -a+constrained_gain(kappas[i],p,kmean,kscale)*G[i] :
                 -a-constrained_gain(kappas[i],p,kmean,kscale)+constrained_gain(kappas[i],p,kmean,kscale)*G[i] for i in eachindex(valid)]
    metric(ids)=begin
        rmse=sqrt(mean((lambdas[ids].-prediction[ids]).^2)); den=sum((lambdas[ids].-mean(lambdas[ids])).^2)
        (;rmse,r2=den>0 ? 1-sum((lambdas[ids].-prediction[ids]).^2)/den : NaN)
    end
    return (;form,a,params=p,kmean,kscale,kappas,lambdas,G,predicted=prediction,
             train=metric(train),test=metric(test),train_indices=train,test_indices=test)
end

constrained_effective_gain(kappa,fit)=constrained_gain(kappa,fit.params,fit.kmean,fit.kscale)

function compare_constrained_couplings(scan)
    fits=[constrained_susceptibility_fit(scan,:additive),constrained_susceptibility_fit(scan,:difference)]
    selected=fits[argmin([f.test.rmse for f in fits])]
    return (;fits,selected)
end

function infer_delayed_coefficients(obs,coupling,kappa,susceptibility)
    a=susceptibility.a; Lambda=constrained_effective_gain(kappa,susceptibility)
    G0=real(delay_transform(coupling,obs.lambda0)); mu=susceptibility.form==:additive ? Lambda-a : -a
    if !(isfinite(obs.dchar)&&isfinite(obs.lambdachar)&&obs.dchar>0)
        return (;form=susceptibility.form,a,Lambda,mu,b=NaN,c=NaN,G0,Gchar=NaN,physical=false,dunstable=NaN,dsn=NaN)
    end
    z=obs.dchar^2; Gc=real(delay_transform(coupling,obs.lambdachar))
    Q=susceptibility.form==:additive ? obs.lambdachar+a-Lambda*Gc : obs.lambdachar+a+Lambda-Lambda*Gc
    v=-(Q+3mu)/2                 # v = c*z^2
    u=v-mu                       # u = b*z
    b=u/z; c=v/z^2
    disc=b^2+4c*mu; zm=c>0&&disc>=0 ? (b-sqrt(disc))/(2c) : NaN
    return (;form=susceptibility.form,a,Lambda,mu,b,c,G0,Gchar=Gc,physical=(a>0&&Lambda>=0&&b>0&&c>0&&mu<0&&disc>=0),
             dunstable=isfinite(zm)&&zm>0 ? sqrt(zm) : NaN,
             dsn=b>0&&c>0 ? sqrt(b/(2c)) : NaN)
end

function periodic_coefficients(mu,b,c,period)
    omega=2pi/period
    target=[mu,b/omega^2,-c/omega^4]
    M=[-1.0 2.0 -3.0; 1/6 -4/3 9/2; -1/120 4/15 -81/40]
    A,B,C=M\target
    return (;A,B,C,omega)
end

periodic_drift(psi,p)= -p.A*sin(psi)+p.B*sin(2psi)-p.C*sin(3psi)
periodic_potential(psi,p)= -p.A*cos(psi)+(p.B/2)*cos(2psi)-(p.C/3)*cos(3psi)

function periodic_equilibria(p;n=4001)
    xs=collect(range(-pi,pi,length=n)); f=periodic_drift.(xs,Ref(p)); roots=Float64[]
    for i in 1:n-1
        if f[i]==0 || f[i]*f[i+1]<0
            lo,hi=xs[i],xs[i+1]
            for _ in 1:50
                mid=(lo+hi)/2
                periodic_drift(lo,p)*periodic_drift(mid,p)<=0 ? (hi=mid) : (lo=mid)
            end
            r=(lo+hi)/2; all(abs(wrap_phase(r-q))>1e-4 for q in roots)&&push!(roots,r)
        end
    end
    rows=NamedTuple[]
    for r in roots
        h=1e-5; deriv=(periodic_drift(r+h,p)-periodic_drift(r-h,p))/(2h)
        push!(rows,(;phase=r,delay=r/p.omega,eigenvalue=deriv,stable=deriv<0))
    end
    return rows
end

function graph_mode_transform(coupling,m,lambda)
    N=length(coupling.neighbors); acc=0.0+0im
    for i in 1:N,q in eachindex(coupling.neighbors[i])
        j=coupling.neighbors[i][q]; delta=j-i
        # unwrap to the shortest signed displacement on the ring
        delta> N÷2 && (delta-=N); delta< -N÷2 && (delta+=N)
        acc += coupling.weights[i][q]*exp(-2pi*im*m*delta/N)*exp(-lambda*coupling.taus[i][q])
    end
    return acc/N
end

function delayed_mode_root(fprime,Lambda,coupling,m;difference=true,initial=fprime)
    lambda=ComplexF64(initial)
    for _ in 1:80
        H=graph_mode_transform(coupling,m,lambda)
        F=lambda-fprime+(difference ? Lambda : 0)-Lambda*H
        epsc=1e-6
        deriv=1-Lambda*(graph_mode_transform(coupling,m,lambda+epsc)-H)/epsc
        abs(deriv)<1e-12&&break
        step=F/deriv; lambda-=step
        abs(step)<1e-10&&break
    end
    return lambda
end

function graph_spectrum(coeff,coupling,N)
    rows=NamedTuple[]
    fp=coeff.mu
    for m in 0:fld(N,2)
        lam=delayed_mode_root(fp,coeff.Lambda,coupling,m;difference=coeff.form==:difference,initial=fp)
        push!(rows,(;mode=m,growth_rate=real(lam),frequency=imag(lam)))
    end
    return rows
end

function reduced_graph_nucleus(coeff,coupling,N,width_fraction,amplitude,rcfg)
    d=zeros(Float64,N); width=max(1,round(Int,width_fraction*N)); lo=fld(N-width,2)+1
    d[lo:lo+width-1].=amplitude
    for _ in 1:rcfg.nucleus_steps
        rhs=coeff.mu.*d .+ coeff.b.*d.^3 .- coeff.c.*d.^5
        for i in 1:N
            avg=sum(coupling.weights[i][q]*d[coupling.neighbors[i][q]] for q in eachindex(coupling.neighbors[i]))
            rhs[i]+=coeff.Lambda*(avg-d[i])
        end
        d .+= rcfg.nucleus_dt.*rhs
        all(isfinite,d)||break
    end
    threshold=isfinite(coeff.dunstable) ? coeff.dunstable : amplitude/2
    return (;final_character_fraction=mean(abs.(d).>threshold),final_mean_abs=mean(abs.(d)),final_max=maximum(abs.(d)))
end

function nucleus_phase_diagram(coeff,coupling,N,rcfg)
    rows=NamedTuple[]; base=isfinite(coeff.dunstable) ? coeff.dunstable : 1.0
    for width in rcfg.nucleus_widths,factor in rcfg.nucleus_amplitude_factors
        amplitude=factor*base; out=reduced_graph_nucleus(coeff,coupling,N,width,amplitude,rcfg)
        push!(rows,(;width_fraction=width,amplitude_factor=factor,amplitude,
                    final_character_fraction=out.final_character_fraction,
                    final_mean_abs=out.final_mean_abs,final_max=out.final_max))
    end
    return rows
end

function plot_revised_analytical(comparison,coeff,periodic,eqs,spectrum,nucleus,outpath)
    sus=comparison.selected
    p1=scatter(sus.kappas,sus.lambdas;label="measured",xlabel="κ",ylabel="λ₀",title="Constrained delayed susceptibility")
    for f in comparison.fits
        plot!(p1,f.kappas,f.predicted;label="$(f.form), test RMSE=$(round(f.test.rmse,digits=3))",lw=f.form==sus.form ? 3 : 1.5)
    end
    phases=range(-pi,pi,length=500); drift=periodic_drift.(phases,Ref(periodic))
    p2=plot(phases,drift;label="periodic drift",xlabel="ψ",ylabel="dψ/dt",title="Periodic tri-stable normal form")
    hline!(p2,[0.0];ls=:dash,label=false)
    for e in eqs; scatter!(p2,[e.phase],[0.0];marker=e.stable ? :circle : :xcross,label=false); end
    p3=plot(getproperty.(spectrum,:mode),getproperty.(spectrum,:growth_rate);marker=:circle,label="Re λ",xlabel="ring mode",ylabel="growth rate",title="Exact nonlocal delayed spectrum")
    widths=sort(unique(getproperty.(nucleus,:width_fraction))); factors=sort(unique(getproperty.(nucleus,:amplitude_factor)))
    Z=fill(NaN,length(factors),length(widths)); lookup=Dict((r.amplitude_factor,r.width_fraction)=>r.final_character_fraction for r in nucleus)
    for (i,f) in enumerate(factors),(j,w) in enumerate(widths); Z[i,j]=lookup[(f,w)]; end
    p4=heatmap(widths,factors,Z;xlabel="perturbed ring fraction",ylabel="amplitude / dᵤ",colorbar_title="final character fraction",title="Critical-nucleus diagram")
    plt=plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm,bottom_margin=12Plots.mm)
    savefig(plt,outpath); return outpath
end

function run_revised_analytical_validation(;cfg=CFG,case=analytical_validation_case(),
        kappas=QUICK ? [case.kappa-.005,case.kappa,case.kappa+.005] : collect(case.kappa-.01:.0025:case.kappa+.01),
        outdir=joinpath(FIGURES_DIR,"revised_analytical_validation"),rcfg=RevisedAnalyticalConfig())
    mkpath(outdir); lc=build_limit_cycle(cfg); baseline=simulate_case(cfg,case,lc;keep_full=true)
    obs=delay_relaxation_observables(baseline)
    scan=NamedTuple[]
    for kappa in kappas
        ccase=merge(case,(;kappa=Float64(kappa))); result=simulate_case(cfg,ccase,lc;keep_full=true)
        o=delay_relaxation_observables(result); G0=isfinite(o.lambda0) ? real(delay_transform(result.coupling,o.lambda0)) : NaN
        push!(scan,(;kappa=Float64(kappa),lambda0=o.lambda0,G0,dchar=o.dchar,lambdachar=o.lambdachar))
    end
    # Retain at least one residual degree of freedom; the quick three-point scan
    # therefore fits only a linear gain instead of an exact quadratic interpolant.
    comparison=compare_constrained_couplings(scan); susceptibility=comparison.selected
    coeff=infer_delayed_coefficients(obs,baseline.coupling,case.kappa,susceptibility)
    periodic=periodic_coefficients(coeff.mu,coeff.b,coeff.c,median(diff(baseline.chorus)))
    eqs=periodic_equilibria(periodic); spectrum=graph_spectrum(coeff,baseline.coupling,cfg.N)
    nucleus=coeff.physical ? nucleus_phase_diagram(coeff,baseline.coupling,cfg.N,rcfg) : NamedTuple[]
    write_namedtuples_csv(joinpath(outdir,"susceptibility_scan.csv"),scan)
    write_namedtuples_csv(joinpath(outdir,"periodic_equilibria.csv"),eqs)
    write_namedtuples_csv(joinpath(outdir,"graph_spectrum.csv"),spectrum)
    write_namedtuples_csv(joinpath(outdir,"critical_nucleus.csv"),nucleus)
    fitrows=[(;form=f.form,a=f.a,test_rmse=f.test.rmse,test_r2=f.test.r2,
               train_rmse=f.train.rmse,train_r2=f.train.r2,selected=f.form==susceptibility.form) for f in comparison.fits]
    write_namedtuples_csv(joinpath(outdir,"coupling_model_comparison.csv"),fitrows)
    summary=[(;form=coeff.form,a=coeff.a,Lambda=coeff.Lambda,mu=coeff.mu,b=coeff.b,c=coeff.c,
              physical=coeff.physical,d_unstable=coeff.dunstable,d_fold=coeff.dsn,
              susceptibility_test_rmse=susceptibility.test.rmse,susceptibility_test_r2=susceptibility.test.r2,
              A=periodic.A,B=periodic.B,C=periodic.C)]
    write_namedtuples_csv(joinpath(outdir,"revised_parameters.csv"),summary)
    !isempty(nucleus)&&plot_revised_analytical(comparison,coeff,periodic,eqs,spectrum,nucleus,joinpath(outdir,"revised_analytical_validation.png"))
    return (;baseline,obs,comparison,susceptibility,coeff,periodic,eqs,spectrum,nucleus,outdir)
end
