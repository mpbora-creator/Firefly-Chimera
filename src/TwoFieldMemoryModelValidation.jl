module TwoFieldMemoryModelValidation

using LinearAlgebra, Statistics, Plots
const SPC = Main.SarfatiPelegChimera

Base.@kwdef struct TwoFieldConfig
    beta_values::Vector{Float64} = [0.25,0.5,1.0,2.0,4.0,8.0]
    gamma_values::Vector{Float64} = [0.01,0.05,0.1,0.25,0.5,1.0,2.0]
    memory_period_factors::Vector{Float64} = [0.5,1.0,2.0,4.0]
    train_fraction::Float64 = 0.70
    ridge::Float64 = 1e-8
    min_rows::Int = 40
    max_mode::Int = SPC.QUICK ? 8 : 50
    root_real_seeds::Int = SPC.QUICK ? 5 : 9
    root_imag_seeds::Int = SPC.QUICK ? 17 : 41
    root_imag_max::Float64 = 12.0
    enforce_tristability::Bool = true
end

memory_function(d,beta,gamma)=d*(1+beta*d^2)/(1+gamma*d^4)

function memory_derivative(d,beta,gamma)
    h=1e-6*max(1.0,abs(d))
    return (memory_function(d+h,beta,gamma)-memory_function(d-h,beta,gamma))/(2h)
end

function reconstruct_memory(D,chorus,beta,gamma,tau_m)
    N,C=size(D); M=zeros(Float64,N,C)
    for i in 1:N
        first=findfirst(isfinite,view(D,i,:)); first===nothing&&continue
        M[i,first]=memory_function(D[i,first],beta,gamma)
        for n in first:C-1
            dt=chorus[n+1]-chorus[n]; alpha=1-exp(-dt/tau_m)
            d=isfinite(D[i,n]) ? D[i,n] : (n>1 ? D[i,n-1] : 0.0)
            M[i,n+1]=(1-alpha)*M[i,n]+alpha*memory_function(d,beta,gamma)
        end
    end
    return M
end

function two_field_design(result,beta,gamma,tau_m)
    D=result.delay_tracks.delay_matrix; chorus=result.chorus; coupling=result.coupling
    N,C=size(D); M=reconstruct_memory(D,chorus,beta,gamma,tau_m)
    X=NTuple{4,Float64}[]; y=Float64[]; cycles=Int[]; oscillator=Int[]
    for n in 2:C-1,i in 1:N
        d0,d1=Float64(D[i,n]),Float64(D[i,n+1]); dt=chorus[n+1]-chorus[n]
        isfinite(d0)&&isfinite(d1)&&dt>0||continue
        dm,dp=Float64(D[mod1(i-1,N),n]),Float64(D[mod1(i+1,N),n])
        isfinite(dm)&&isfinite(dp)||continue
        graph=0.0; valid=true
        for q in eachindex(coupling.neighbors[i])
            j=coupling.neighbors[i][q]
            dj=SPC.interp_delay_track(chorus[n]-coupling.taus[i][q],chorus,view(D,j,:))
            if !isfinite(dj); valid=false; break; end
            graph+=coupling.weights[i][q]*(dj-d0)
        end
        valid||continue
        push!(X,(-d0,M[i,n],dm-2d0+dp,graph)); push!(y,(d1-d0)/dt)
        push!(cycles,n); push!(oscillator,i)
    end
    Xm=isempty(X) ? zeros(0,4) : [X[i][j] for i in eachindex(X),j in 1:4]
    return (;X=Xm,y,cycles,oscillator,memory=M,D,chorus)
end

function nonnegative_subset_fit(X,y,ridge)
    p=size(X,2); best=fill(0.0,p); bestloss=Inf
    for mask in 1:(2^p-1)
        cols=[j for j in 1:p if mask&(1<<(j-1))!=0]
        beta=(X[:,cols]'X[:,cols]+ridge*Matrix{Float64}(I,length(cols),length(cols)))\(X[:,cols]'y)
        all(beta.>=0)||continue
        full=zeros(p); full[cols].=beta; loss=mean((y.-X*full).^2)
        if loss<bestloss; bestloss=loss; best=full; end
    end
    return best
end

function nonnegative_fit_required_first(X,y,ridge)
    p=size(X,2); best=nothing; bestloss=Inf
    # The first coefficient is chi and must be active; optional spatial terms
    # are selected by subset enumeration.
    for mask in 0:(2^(p-1)-1)
        cols=vcat(1,[j for j in 2:p if mask&(1<<(j-2))!=0])
        beta=(X[:,cols]'X[:,cols]+ridge*Matrix{Float64}(I,length(cols),length(cols)))\(X[:,cols]'y)
        all(beta.>=0)&&beta[1]>1e-10||continue
        full=zeros(p); full[cols].=beta; loss=mean((y.-X*full).^2)
        if loss<bestloss; bestloss=loss; best=full; end
    end
    return best
end

function metrics(y,pred)
    rmse=sqrt(mean((y.-pred).^2)); den=sum((y.-mean(y)).^2)
    return (;rmse,r2=den>0 ? 1-sum((y.-pred).^2)/den : NaN,
             correlation=std(y)>0&&std(pred)>0 ? cor(y,pred) : NaN)
end

function fit_candidate(result,beta,gamma,tau_m,cfg)
    data=two_field_design(result,beta,gamma,tau_m); n=length(data.y)
    n>=cfg.min_rows||return nothing
    order=sortperm(data.cycles); cut=clamp(round(Int,cfg.train_fraction*n),1,n-1)
    train,test=order[1:cut],order[cut+1:end]
    coef=nonnegative_subset_fit(data.X[train,:],data.y[train],cfg.ridge)
    return (;beta,gamma,tau_m,coefficients=coef,data,train,test,
             train_metrics=metrics(data.y[train],data.X[train,:]*coef),
             test_metrics=metrics(data.y[test],data.X[test,:]*coef))
end

function measured_character_delay(result)
    vals=abs.(result.delay_tracks.stable_delay_mean)
    !isempty(vals)&&return median(vals)
    vals=abs.(result.classif.character_delays)
    return isempty(vals) ? NaN : median(vals)
end

function fit_tristable_candidate(result,beta,gamma,tau_m,dtarget,cfg)
    isfinite(dtarget)&&dtarget>0||return nothing
    data=two_field_design(result,beta,gamma,tau_m); n=length(data.y)
    n>=cfg.min_rows||return nothing
    z=dtarget^2; ratio=(1+beta*z)/(1+gamma*z^2)
    gp=memory_derivative(dtarget,beta,gamma)
    # a=chi*ratio gives the measured nonzero equilibrium exactly. The two
    # inequalities guarantee stable chorus and stable character branches.
    ratio>1&&ratio>gp||return nothing
    X=hcat(data.X[:,2].-ratio.*(-data.X[:,1]),data.X[:,3],data.X[:,4])
    order=sortperm(data.cycles); cut=clamp(round(Int,cfg.train_fraction*n),1,n-1)
    train,test=order[1:cut],order[cut+1:end]
    reduced=nonnegative_fit_required_first(X[train,:],data.y[train],cfg.ridge); reduced===nothing&&return nothing
    chi,Dd,Lambda=reduced; a=ratio*chi; coef=[a,chi,Dd,Lambda]
    return (;beta,gamma,tau_m,coefficients=coef,data,train,test,target_delay=dtarget,
             equilibrium_ratio=ratio,character_slope=gp,
             train_metrics=metrics(data.y[train],X[train,:]*reduced),
             test_metrics=metrics(data.y[test],X[test,:]*reduced),joint_X=X)
end

function calibrate_two_field(result;cfg=TwoFieldConfig())
    T=median(diff(result.chorus)); fits=NamedTuple[]; dtarget=measured_character_delay(result)
    for beta in cfg.beta_values,gamma in cfg.gamma_values,factor in cfg.memory_period_factors
        fit=cfg.enforce_tristability ? fit_tristable_candidate(result,beta,gamma,factor*T,dtarget,cfg) :
                                      fit_candidate(result,beta,gamma,factor*T,cfg)
        fit===nothing||push!(fits,fit)
    end
    isempty(fits)&&error("no two-field candidate had enough valid observations")
    selected=fits[argmin([f.test_metrics.rmse for f in fits])]
    return (;fits,selected,chorus_period=T,target_delay=dtarget)
end

function homogeneous_equilibria(fit)
    a,chi=fit.coefficients[1],fit.coefficients[2]; beta,gamma=fit.beta,fit.gamma
    rows=[(;delay=0.0,memory=0.0,stable=a>chi,eigen_trace=-a-1/fit.tau_m,
            eigen_det=(a-chi)/fit.tau_m)]
    A=a*gamma; B=-chi*beta; C=a-chi; disc=B^2-4A*C
    if A!=0&&disc>=0
        for z in ((-B+sqrt(disc))/(2A),(-B-sqrt(disc))/(2A))
            z>0||continue
            for d in (-sqrt(z),sqrt(z))
                gp=memory_derivative(d,beta,gamma); det=(a-chi*gp)/fit.tau_m
                push!(rows,(;delay=d,memory=memory_function(d,beta,gamma),stable=det>0,
                            eigen_trace=-a-1/fit.tau_m,eigen_det=det))
            end
        end
    end
    sort!(rows,by=x->x.delay); return rows
end

function transform(coupling,mode,lambda)
    N=length(coupling.neighbors); value=0.0+0im
    for i in 1:N,q in eachindex(coupling.neighbors[i])
        j=coupling.neighbors[i][q]; delta=j-i
        delta>N÷2&&(delta-=N); delta< -N÷2&&(delta+=N)
        value+=coupling.weights[i][q]*exp(-2pi*im*mode*delta/N-lambda*coupling.taus[i][q])
    end
    value/N
end

function characteristic(lambda,eq,fit,coupling,mode,N)
    a,chi,Dd,Lambda=fit.coefficients; k=2pi*mode/N; Dm=0.0
    gp=memory_derivative(eq.delay,fit.beta,fit.gamma)
    return (lambda+a+Dd*k^2+Lambda*(1-transform(coupling,mode,lambda)))*(fit.tau_m*lambda+1+Dm*k^2)-chi*gp
end

function characteristic_root(seed,eq,fit,coupling,mode,N)
    z=ComplexF64(seed)
    for _ in 1:100
        F=characteristic(z,eq,fit,coupling,mode,N); h=1e-6
        der=(characteristic(z+h,eq,fit,coupling,mode,N)-F)/h
        abs(der)<1e-12&&return nothing
        step=F/der; z-=step; abs(step)<1e-9&&break
    end
    abs(characteristic(z,eq,fit,coupling,mode,N))<1e-6 ? z : nothing
end

function two_field_spectrum(equilibria,fit,coupling,N,cfg)
    rows=NamedTuple[]
    for eq in equilibria,m in 0:min(cfg.max_mode,N÷2)
        roots=ComplexF64[]
        for x in range(-2,.5,length=cfg.root_real_seeds),y in range(-cfg.root_imag_max,cfg.root_imag_max,length=cfg.root_imag_seeds)
            r=characteristic_root(x+im*y,eq,fit,coupling,m,N); r===nothing&&continue
            all(abs(r-q)>1e-5 for q in roots)&&push!(roots,r)
        end
        sort!(roots,by=real,rev=true)
        for (branch,r) in enumerate(roots)
            push!(rows,(;equilibrium_delay=eq.delay,mode=m,branch,growth_rate=real(r),frequency=imag(r),dominant=branch==1))
        end
    end
    rows
end

function write_csv(path,rows)
    isempty(rows)&&return
    open(path,"w") do io
        names=propertynames(first(rows)); println(io,join(string.(names),","))
        for r in rows; println(io,join([string(getproperty(r,n)) for n in names],",")); end
    end
end

function run(;case=SPC.analytical_validation_case(),cfg=TwoFieldConfig(),
        outdir=joinpath(SPC.FIGURES_DIR,"two_field_memory_validation"))
    mkpath(outdir); lc=SPC.build_limit_cycle(SPC.CFG)
    result=SPC.simulate_case(SPC.CFG,case,lc;keep_full=true); calibration=calibrate_two_field(result;cfg)
    fit=calibration.selected; equilibria=homogeneous_equilibria(fit)
    spectrum=two_field_spectrum(equilibria,fit,result.coupling,SPC.CFG.N,cfg)
    coeff=fit.coefficients
    write_csv(joinpath(outdir,"candidate_models.csv"),[(;beta=f.beta,gamma=f.gamma,tau_m=f.tau_m,
        target_delay=hasproperty(f,:target_delay) ? f.target_delay : NaN,
        test_rmse=f.test_metrics.rmse,test_r2=f.test_metrics.r2) for f in calibration.fits])
    write_csv(joinpath(outdir,"selected_parameters.csv"),[(;a=coeff[1],chi=coeff[2],D_d=coeff[3],Lambda=coeff[4],
        beta=fit.beta,gamma=fit.gamma,tau_m=fit.tau_m,test_rmse=fit.test_metrics.rmse,
        test_r2=fit.test_metrics.r2,test_correlation=fit.test_metrics.correlation,
        target_character_delay=calibration.target_delay,tristability_enforced=cfg.enforce_tristability)])
    write_csv(joinpath(outdir,"homogeneous_equilibria.csv"),equilibria)
    write_csv(joinpath(outdir,"two_field_spectrum.csv"),spectrum)
    pred=(hasproperty(fit,:joint_X) ? fit.joint_X[fit.test,:]*[coeff[2],coeff[3],coeff[4]] : fit.data.X[fit.test,:]*coeff); obs=fit.data.y[fit.test]
    p1=scatter(obs,pred;ms=2,alpha=.4,label=false,xlabel="measured d-dot",ylabel="predicted d-dot",title="Two-field held-out prediction")
    lim=collect(extrema(vcat(obs,pred))); plot!(p1,lim,lim;ls=:dash,label="identity")
    p2=scatter(vec(fit.data.D),vec(fit.data.memory);ms=1,alpha=.15,label=false,xlabel="delay d",ylabel="memory m",title="Reconstructed slow memory")
    p3=scatter(getproperty.(equilibria,:delay),zeros(length(equilibria));marker_z=Float64.(getproperty.(equilibria,:stable)),label=false,
        xlabel="equilibrium delay",ylabel="",yticks=false,title="Predicted equilibria")
    dom=filter(x->x.dominant,spectrum); p4=plot(;xlabel="mode",ylabel="dominant growth",title="Two-field distributed-delay spectrum")
    for d in unique(getproperty.(dom,:equilibrium_delay)); r=filter(x->x.equilibrium_delay==d,dom); plot!(p4,getproperty.(r,:mode),getproperty.(r,:growth_rate);marker=:circle,label="d=$(round(d,digits=2))"); end
    savefig(plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm),joinpath(outdir,"two_field_memory_validation.png"))
    return (;result,calibration,fit,equilibria,spectrum,outdir)
end

end # module
