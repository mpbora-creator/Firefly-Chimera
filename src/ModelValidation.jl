# Numerical validation tools for the reduced delay-field model in Section 5.
# These routines deliberately use the event-defined d_i(n) field produced by
# Diagnostics.jl and the same neighbours, weights, and delays as the simulator.

Base.@kwdef struct ValidationConfig
    ridge::Float64 = 1.0
    train_fraction::Float64 = 0.70
    min_rows::Int = 30
    branch_block_fraction::Float64 = 0.15
    branch_phase_offset::Float64 = 0.70
    regime_delay_threshold::Float64 = CHARACTER_TRACK_MIN_DELAY
    min_regime_rows::Int = 12
end

const REDUCED_FEATURES = (:constant, :d, :d2, :d3, :d4, :d5,
                          :d_previous, :graph_difference, :sin1, :sin2)
const CIRCULAR_FEATURES = (:constant, :cos_phase, :sin_phase,
                           :cos_previous, :sin_previous,
                           :graph_sin1, :graph_sin2, :current_character)

wrap_phase(x) = mod(x + pi, 2pi) - pi

function interp_delay_track(t::Real, chorus, row)
    ids = findall(isfinite, row)
    isempty(ids) && return NaN
    tv = chorus[ids]
    k = searchsortedfirst(tv, t)
    if k <= 1
        return abs(t-tv[1]) <= (length(tv) > 1 ? tv[2]-tv[1] : Inf) ? Float64(row[ids[1]]) : NaN
    elseif k > length(tv)
        return abs(t-tv[end]) <= (length(tv) > 1 ? tv[end]-tv[end-1] : Inf) ? Float64(row[ids[end]]) : NaN
    end
    i0, i1 = ids[k-1], ids[k]
    t0, t1 = chorus[i0], chorus[i1]
    t1 == t0 && return Float64(row[i0])
    q = (t-t0)/(t1-t0)
    return (1-q)*Float64(row[i0]) + q*Float64(row[i1])
end

function validation_design(result)
    D = result.delay_tracks.delay_matrix
    chorus = result.chorus
    coupling = result.coupling
    N, C = size(D)
    X = Vector{NTuple{10,Float64}}()
    y = Float64[]
    oscillator = Int[]
    cycle = Int[]
    C < 3 && return (; X=zeros(0,10), y, oscillator, cycle)
    omega = 2pi/median(diff(chorus))

    for n in 2:C-1, i in 1:N
        dprev, d0, d1 = Float64(D[i,n-1]), Float64(D[i,n]), Float64(D[i,n+1])
        isfinite(dprev) && isfinite(d0) && isfinite(d1) || continue
        graphdiff = 0.0
        sin1 = 0.0
        sin2 = 0.0
        valid = true
        for q in eachindex(coupling.neighbors[i])
            j = coupling.neighbors[i][q]
            dj = interp_delay_track(chorus[n]-coupling.taus[i][q], chorus, view(D,j,:))
            if !isfinite(dj)
                valid = false
                break
            end
            delta = dj-d0
            w = coupling.weights[i][q]
            graphdiff += w*delta
            sin1 += w*sin(omega*delta)
            sin2 += w*sin(2omega*delta)
        end
        valid || continue
        push!(X, (1.0,d0,d0^2,d0^3,d0^4,d0^5,dprev,graphdiff,sin1,sin2))
        push!(y, d1)
        push!(oscillator, i); push!(cycle, n)
    end
    Xm = isempty(X) ? zeros(0,10) : [X[i][j] for i in eachindex(X), j in 1:10]
    return (; X=Xm, y, oscillator, cycle, omega)
end

function regression_metrics(y, yhat)
    isempty(y) && return (; rmse=NaN, r2=NaN, correlation=NaN)
    rmse = sqrt(mean((y.-yhat).^2))
    den = sum((y.-mean(y)).^2)
    r2 = den > 0 ? 1-sum((y.-yhat).^2)/den : NaN
    correlation = length(y) > 1 && std(y) > 0 && std(yhat) > 0 ? cor(y,yhat) : NaN
    return (; rmse, r2, correlation)
end

function ridge_fit(X, y, ridge)
    p = size(X,2)
    scales = [std(view(X,:,j)) for j in 1:p]
    scales .= map(x -> isfinite(x) && x > eps() ? x : 1.0, scales)
    Z = X ./ permutedims(scales)
    beta_scaled = (Z'Z + ridge*Matrix{Float64}(I,p,p)) \ (Z'y)
    return beta_scaled ./ scales, beta_scaled
end

function identify_reduced_model(result; cfg=ValidationConfig())
    data = validation_design(result)
    m = length(data.y)
    m >= cfg.min_rows || error("only $m valid delay-field rows; need at least $(cfg.min_rows)")
    order = sortperm(data.cycle)
    cut = clamp(round(Int,cfg.train_fraction*m), 1, m-1)
    train, test = order[1:cut], order[cut+1:end]
    models = [
        (:persistence, [1,2,7]),
        (:symmetric_pitchfork, [2,4,6]),
        (:asymmetric_local, collect(1:7)),
        (:linear_graph_map, collect(1:8)),
        (:harmonic_graph_map, [1,2,3,4,5,6,7,9,10]),
        (:full_event_map, collect(1:10)),
    ]
    fits = NamedTuple[]
    for (name, cols) in models
        beta, beta_scaled = ridge_fit(data.X[train,cols], data.y[train], cfg.ridge)
        pred_train = data.X[train,cols]*beta
        pred_test = data.X[test,cols]*beta
        fullbeta = fill(NaN, length(REDUCED_FEATURES)); fullbeta[cols] .= beta
        fullscaled = fill(NaN, length(REDUCED_FEATURES)); fullscaled[cols] .= beta_scaled
        push!(fits, (; name, columns=Tuple(REDUCED_FEATURES[c] for c in cols), coefficients=fullbeta,
                     scaled_coefficients=fullscaled,
                     train=regression_metrics(data.y[train],pred_train),
                     test=regression_metrics(data.y[test],pred_test)))
    end
    selected = fits[argmin([x.test.rmse for x in fits])]
    return (; data, train, test, fits, selected, full_event=last(fits))
end

function circular_validation_design(result; cfg=ValidationConfig())
    D=result.delay_tracks.delay_matrix; chorus=result.chorus; coupling=result.coupling
    N,C=size(D); C < 3 && return (;X=zeros(0,8),psi=Float64[],regime=Bool[],next_regime=Bool[],cycle=Int[])
    T=median(diff(chorus)); omega=2pi/T
    rows=Vector{NTuple{8,Float64}}(); target=Float64[]; regime=Bool[]; nextregime=Bool[]; cycles=Int[]
    for n in 2:C-1, i in 1:N
        dp,d0,d1=Float64(D[i,n-1]),Float64(D[i,n]),Float64(D[i,n+1])
        isfinite(dp)&&isfinite(d0)&&isfinite(d1) || continue
        pprev=wrap_phase(omega*dp); p0=wrap_phase(omega*d0)
        gs1=0.0; gs2=0.0; valid=true
        for q in eachindex(coupling.neighbors[i])
            j=coupling.neighbors[i][q]
            dj=interp_delay_track(chorus[n]-coupling.taus[i][q],chorus,view(D,j,:))
            if !isfinite(dj); valid=false; break; end
            delta=wrap_phase(omega*(dj-d0)); w=coupling.weights[i][q]
            gs1 += w*sin(delta); gs2 += w*sin(2delta)
        end
        valid || continue
        current_character=abs(d0)>cfg.regime_delay_threshold
        push!(rows,(1.0,cos(p0),sin(p0),cos(pprev),sin(pprev),gs1,gs2,current_character ? 1.0 : 0.0))
        push!(target,wrap_phase(omega*d1))
        push!(regime,current_character); push!(nextregime,abs(d1)>cfg.regime_delay_threshold); push!(cycles,n)
    end
    X=isempty(rows) ? zeros(0,8) : [rows[i][j] for i in eachindex(rows),j in 1:8]
    return (;X,psi=target,regime,next_regime=nextregime,cycle=cycles,omega,T)
end

function circular_metrics(obs,pred)
    err=wrap_phase.(pred.-obs)
    return (; angular_rmse=sqrt(mean(err.^2)), mean_resultant=mean(cos.(err)),
             mean_bias=atan(mean(sin.(err)),mean(cos.(err))))
end

function circular_fit_predict(X,psi,train,test,cols,ridge)
    bc,_=ridge_fit(X[train,cols],cos.(psi[train]),ridge)
    bs,_=ridge_fit(X[train,cols],sin.(psi[train]),ridge)
    pred=atan.(X[test,cols]*bs,X[test,cols]*bc)
    return (;bc,bs,pred)
end

logistic(x)=1/(1+exp(-clamp(x,-35.0,35.0)))

function logistic_ridge_fit(X,y,ridge;iterations=60)
    p=size(X,2); scales=[std(view(X,:,j)) for j in 1:p]
    scales .= map(x->isfinite(x)&&x>eps() ? x : 1.0,scales)
    Z=X ./ permutedims(scales); beta=zeros(p)
    for _ in 1:iterations
        prob=logistic.(Z*beta); w=max.(prob.*(1 .- prob),1e-6)
        grad=Z'*(y.-prob)-ridge*beta
        H=Z'*(Z.*w)+ridge*Matrix{Float64}(I,p,p)
        step=H\grad; beta .+= step
        norm(step)<1e-8 && break
    end
    return beta./scales
end

function binary_metrics(y,p)
    pp=clamp.(p,1e-9,1-1e-9); yy=Float64.(y)
    logloss=-mean(yy.*log.(pp)+(1 .- yy).*log.(1 .- pp))
    brier=mean((pp.-yy).^2)
    pos=findall(y); neg=findall(.!y)
    auc=isempty(pos)||isempty(neg) ? NaN : mean([p[i]>p[j] ? 1.0 : p[i]==p[j] ? 0.5 : 0.0 for i in pos,j in neg])
    hard=pp.>=0.5
    recall=isempty(pos) ? NaN : mean(hard[pos]); specificity=isempty(neg) ? NaN : mean(.!hard[neg])
    return (;logloss,brier,auc,recall,specificity,prevalence=mean(yy))
end

function probabilistic_switching_fit(data,train,test,cfg)
    cols=collect(1:8); ytrain=Float64.(data.next_regime[train])
    beta=logistic_ridge_fit(data.X[train,cols],ytrain,cfg.ridge)
    probability=logistic.(data.X[test,cols]*beta)
    transition=binary_metrics(data.next_regime[test],probability)

    # Conditional delay maps are trained on the role at n+1, not the current role.
    conditional=Dict{Bool,Any}()
    globalfit=circular_fit_predict(data.X,data.psi,train,test,collect(1:7),cfg.ridge)
    predictions=Dict{Bool,Vector{Float64}}()
    for state in (false,true)
        tr=[i for i in train if data.next_regime[i]==state]
        if length(tr)>=cfg.min_regime_rows
            f=circular_fit_predict(data.X,data.psi,tr,test,collect(1:7),cfg.ridge)
            predictions[state]=f.pred; conditional[state]=(;rows=length(tr),bc=f.bc,bs=f.bs)
        else
            predictions[state]=globalfit.pred; conditional[state]=(;rows=length(tr),bc=globalfit.bc,bs=globalfit.bs)
        end
    end
    pc=predictions[true]; ph=predictions[false]
    mixcos=(1 .- probability).*cos.(ph)+probability.*cos.(pc)
    mixsin=(1 .- probability).*sin.(ph)+probability.*sin.(pc)
    prediction=atan.(mixsin,mixcos)
    return (;name=:probabilistic_switching_graph,probability,prediction,
             metrics=circular_metrics(data.psi[test],prediction),transition,beta,conditional,test_indices=test)
end

function identify_circular_model(result; cfg=ValidationConfig())
    data=circular_validation_design(result;cfg)
    m=length(data.psi); m>=cfg.min_rows || error("only $m valid circular-map rows")
    order=sortperm(data.cycle); cut=clamp(round(Int,cfg.train_fraction*m),1,m-1)
    train,test=order[1:cut],order[cut+1:end]
    specifications=[(:circular_persistence,collect(1:5)),(:circular_graph,collect(1:7))]
    fits=NamedTuple[]
    for (name,cols) in specifications
        f=circular_fit_predict(data.X,data.psi,train,test,cols,cfg.ridge)
        push!(fits,(;name,columns=Tuple(CIRCULAR_FEATURES[c] for c in cols),
                     metrics=circular_metrics(data.psi[test],f.pred),prediction=f.pred,
                     cosine_coefficients=f.bc,sine_coefficients=f.bs,test_indices=test))
    end

    # Switching model: fit separate circular graph maps to chorus and character rows.
    cols=collect(1:7); pred=fill(NaN,length(test)); regime_fits=Dict{Bool,Any}()
    globalfit=circular_fit_predict(data.X,data.psi,train,test,cols,cfg.ridge)
    for state in (false,true)
        tr=[i for i in train if data.regime[i]==state]
        tepos=[k for (k,i) in enumerate(test) if data.regime[i]==state]
        if length(tr)>=cfg.min_regime_rows && !isempty(tepos)
            te=test[tepos]; f=circular_fit_predict(data.X,data.psi,tr,te,cols,cfg.ridge)
            pred[tepos].=f.pred; regime_fits[state]=(;rows=length(tr),bc=f.bc,bs=f.bs)
        else
            pred[tepos].=globalfit.pred[tepos]
            regime_fits[state]=(;rows=length(tr),bc=globalfit.bc,bs=globalfit.bs)
        end
    end
    push!(fits,(;name=:switching_circular_graph,columns=Tuple(CIRCULAR_FEATURES),
                 metrics=circular_metrics(data.psi[test],pred),prediction=pred,
                 cosine_coefficients=regime_fits,sine_coefficients=regime_fits,test_indices=test))
    probabilistic=probabilistic_switching_fit(data,train,test,cfg)
    push!(fits,(;name=probabilistic.name,columns=Tuple(CIRCULAR_FEATURES),metrics=probabilistic.metrics,
                 prediction=probabilistic.prediction,cosine_coefficients=probabilistic.conditional,
                 sine_coefficients=probabilistic.conditional,test_indices=test))
    selected=fits[argmin([x.metrics.angular_rmse for x in fits])]
    return (;data,train,test,fits,selected,probabilistic)
end

function controlled_phase_profile(N, kind::Symbol, fraction, offset)
    phases = zeros(Float64,N)
    kind == :chorus && return phases
    kind in (:positive,:negative) || error("unknown branch kind $kind")
    width = max(1,round(Int,fraction*N))
    start = fld(N-width,2)+1
    phases[start:start+width-1] .= kind == :positive ? offset : -offset
    return phases
end

function branch_sweep(cfg::SearchConfig, basecase, lc, kappas;
                      vcfg=ValidationConfig(), kinds=(:chorus,:positive,:negative))
    rows = NamedTuple[]
    for kappa in kappas, kind in kinds
        case = merge(basecase,(; kappa=Float64(kappa)))
        phases = controlled_phase_profile(cfg.N,kind,vcfg.branch_block_fraction,vcfg.branch_phase_offset)
        result = simulate_case(cfg,case,lc;keep_full=true,initial_phases=phases)
        D = result.delay_tracks.delay_matrix
        vals = D[isfinite.(D)]
        tailstart = max(1,size(D,2)-max(3,round(Int,0.25size(D,2)))+1)
        tail = D[:,tailstart:end]; tailvals = tail[isfinite.(tail)]
        push!(rows,(; kappa=Float64(kappa),branch=kind,
                    mean_delay=isempty(tailvals) ? NaN : mean(tailvals),
                    mean_abs_delay=isempty(tailvals) ? NaN : mean(abs.(tailvals)),
                    delay_std=isempty(tailvals) ? NaN : std(tailvals),
                    stable_character_fraction=result.summary.stable_character_fraction,
                    synchrony_fraction=result.summary.synchrony_fraction,
                    score=result.summary.score, observations=length(vals)))
    end
    return rows
end

function write_namedtuples_csv(path, rows)
    isempty(rows) && return path
    open(path,"w") do io
        names = propertynames(first(rows)); println(io,join(string.(names),","))
        for row in rows
            println(io,join((replace(string(getproperty(row,n)),','=>';') for n in names),","))
        end
    end
    return path
end

function plot_identification(fit, outpath)
    f = fit.selected; data=fit.data; test=fit.test
    cols = findall(isfinite,f.coefficients)
    pred = data.X[test,cols]*f.coefficients[cols]
    p1 = scatter(data.y[test],pred;ms=2,alpha=.45,label="samples",xlabel="measured d(n+1)",ylabel="predicted d(n+1)",title="Selected map: $(f.name)")
    lim = collect(extrema(vcat(data.y[test],pred))); plot!(p1,lim,lim;label="identity",lw=2,ls=:dash)
    active=findall(isfinite,f.scaled_coefficients)
    names=[string(REDUCED_FEATURES[i]) for i in active]; vals=f.scaled_coefficients[active]
    p2 = bar(eachindex(vals),vals;xticks=(collect(eachindex(vals)),names),legend=false,xrotation=35,ylabel="scaled coefficient",title="Selected-map coefficients")
    errors=[x.test.rmse for x in fit.fits]; modelnames=string.([x.name for x in fit.fits])
    p3 = bar(eachindex(errors),errors;xticks=(collect(eachindex(errors)),modelnames),legend=false,xrotation=25,ylabel="test RMSE",title="Term ablation")
    residual=data.y[test]-pred
    p4=scatter(data.X[test,2],residual;ms=2,alpha=.4,legend=false,xlabel="d(n)",ylabel="residual",title="Residual structure")
    hline!(p4,[0.0];color=:black,ls=:dash)
    plt=plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm,bottom_margin=14Plots.mm)
    savefig(plt,outpath); return outpath
end

function plot_branch_sweep(rows,outpath)
    p1=plot(;xlabel="κ",ylabel="signed tail mean",title="Controlled initial conditions")
    p2=plot(;xlabel="κ",ylabel="mean |delay|",title="Delay-branch amplitude")
    p3=plot(;xlabel="κ",ylabel="off-chorus flash fraction",title="Character activity")
    for kind in unique(getproperty.(rows,:branch))
        r=sort(filter(x->x.branch==kind,rows),by=x->x.kappa)
        plot!(p1,getproperty.(r,:kappa),getproperty.(r,:mean_delay);marker=:circle,label=string(kind))
        plot!(p2,getproperty.(r,:kappa),getproperty.(r,:mean_abs_delay);marker=:circle,label=string(kind))
        plot!(p3,getproperty.(r,:kappa),1 .- getproperty.(r,:synchrony_fraction);marker=:circle,label=string(kind))
    end
    plt=plot(p1,p2,p3;layout=(1,3),size=(1400,430),margin=7Plots.mm)
    savefig(plt,outpath); return outpath
end

function plot_circular_identification(fit,outpath)
    f=fit.selected; obs=fit.data.psi[fit.test]; pred=f.prediction
    err=wrap_phase.(pred.-obs)
    p1=scatter(obs,pred;ms=2,alpha=.4,label="events",xlabel="measured phase delay",ylabel="predicted phase delay",title="Selected: $(f.name)")
    lo,hi=extrema(vcat(obs,pred)); pad=max(0.05,0.08*(hi-lo)); lim=[max(-pi,lo-pad),min(pi,hi+pad)]
    plot!(p1,lim,lim;ls=:dash,lw=2,label="identity",xlims=Tuple(lim),ylims=Tuple(lim))
    rms=[x.metrics.angular_rmse for x in fit.fits]; labels=string.([x.name for x in fit.fits])
    p2=bar(eachindex(rms),rms;xticks=(collect(eachindex(rms)),labels),xrotation=25,legend=false,ylabel="angular RMSE",title="Circular-model comparison")
    p3=histogram(err;bins=30,normalize=:pdf,legend=false,xlabel="wrapped prediction error",ylabel="density",title="Circular residuals")
    states=fit.data.regime[fit.test]
    roleerr=[any(.!states) ? mean(abs.(err[.!states])) : NaN,
             any(states) ? mean(abs.(err[states])) : NaN]
    p4=bar(1:2,roleerr;xticks=([1,2],["chorus","character"]),legend=false,ylabel="mean |angular error|",title="Error by current role")
    plt=plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm,bottom_margin=14Plots.mm)
    savefig(plt,outpath); return outpath
end

function plot_probabilistic_switching(fit,outpath)
    f=fit.probabilistic; test=fit.test; y=fit.data.next_regime[test]; p=f.probability
    thresholds=collect(0.0:0.02:1.0); pos=max(count(y),1); neg=max(count(.!y),1)
    tpr=[count((p.>=q).&y)/pos for q in thresholds]
    fpr=[count((p.>=q).&(.!y))/neg for q in thresholds]
    p1=plot(fpr,tpr;label="AUC=$(round(f.transition.auc,digits=3))",xlabel="false-positive rate",ylabel="true-positive rate",title="Role-transition ROC",lw=2)
    plot!(p1,[0,1],[0,1];ls=:dash,label="chance")
    edges=collect(0.0:0.1:1.0); centers=0.5.*(edges[1:end-1].+edges[2:end]); observed=fill(NaN,length(centers)); counts=zeros(Int,length(centers))
    for k in eachindex(centers)
        lower=p.>=edges[k]; upper=k==length(centers) ? p.<=edges[k+1] : p.<edges[k+1]
        ids=findall(lower.&upper)
        if !isempty(ids); observed[k]=mean(y[ids]); counts[k]=length(ids); end
    end
    valid=isfinite.(observed); p2=scatter(centers[valid],observed[valid];ms=5,label="bins",xlabel="predicted character probability",ylabel="observed fraction",title="Transition calibration")
    plot!(p2,[0,1],[0,1];ls=:dash,label="ideal")
    p3=histogram(p[.!y];bins=20,alpha=.55,label="next chorus",xlabel="predicted character probability",ylabel="count",title="Probability separation")
    histogram!(p3,p[y];bins=20,alpha=.55,label="next character")
    err=abs.(wrap_phase.(f.prediction.-fit.data.psi[test]))
    p4=scatter(p,err;ms=2,alpha=.4,legend=false,xlabel="predicted character probability",ylabel="|angular error|",title="Uncertainty vs delay error")
    plt=plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=7Plots.mm)
    savefig(plt,outpath); return outpath
end

function run_model_validation(; cfg=CFG, case=SAMPLE_CASE,
        kappas=QUICK ? [case.kappa-0.005,case.kappa,case.kappa+0.005] : collect(case.kappa-0.01:0.0025:case.kappa+0.01),
        outdir=joinpath(FIGURES_DIR,"model_validation"), run_branches=true)
    mkpath(outdir)
    lc=build_limit_cycle(cfg)
    baseline=simulate_case(cfg,case,lc;keep_full=true)
    fit=identify_reduced_model(baseline)
    circular_fit=identify_circular_model(baseline)
    plot_identification(fit,joinpath(outdir,"reduced_model_identification.png"))
    plot_circular_identification(circular_fit,joinpath(outdir,"circular_model_identification.png"))
    plot_probabilistic_switching(circular_fit,joinpath(outdir,"probabilistic_switching.png"))
    coeffrows=[(; model=x.name,test_rmse=x.test.rmse,test_r2=x.test.r2,test_correlation=x.test.correlation,
                 constant=x.coefficients[1],d=x.coefficients[2],d2=x.coefficients[3],d3=x.coefficients[4],
                 d4=x.coefficients[5],d5=x.coefficients[6],d_previous=x.coefficients[7],
                 graph_difference=x.coefficients[8],sin1=x.coefficients[9],sin2=x.coefficients[10]) for x in fit.fits]
    write_namedtuples_csv(joinpath(outdir,"reduced_model_coefficients.csv"),coeffrows)
    circularrows=[(;model=x.name,angular_rmse=x.metrics.angular_rmse,
                    mean_resultant=x.metrics.mean_resultant,mean_bias=x.metrics.mean_bias) for x in circular_fit.fits]
    write_namedtuples_csv(joinpath(outdir,"circular_model_metrics.csv"),circularrows)
    tr=circular_fit.probabilistic.transition
    write_namedtuples_csv(joinpath(outdir,"transition_metrics.csv"),[(;model=:probabilistic_switching_graph,
        logloss=tr.logloss,brier=tr.brier,auc=tr.auc,recall=tr.recall,
        specificity=tr.specificity,character_prevalence=tr.prevalence)])
    branches=NamedTuple[]
    if run_branches
        branches=branch_sweep(cfg,case,lc,kappas)
        write_namedtuples_csv(joinpath(outdir,"branch_sweep.csv"),branches)
        plot_branch_sweep(branches,joinpath(outdir,"branch_sweep.png"))
    end
    return (; baseline,fit,circular_fit,branches,outdir)
end
