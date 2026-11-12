module TwoFieldBasinBoundaryValidation

using Statistics, Random, DelimitedFiles, Plots
const SP = Main.SarfatiPelegChimera

Base.@kwdef struct BasinConfig
    settle_time::Float64 = SP.QUICK ? 90.0 : 260.0
    observe_time::Float64 = SP.QUICK ? 150.0 : 360.0
    block_fractions::Vector{Float64} = SP.QUICK ? [0.30] : [0.15, 0.30, 0.45]
    # Geometric phase is only the actuator. Thresholds are reported using the
    # measured first-flash delay, since relaxation-oscillator phase is nonuniform.
    kick_amplitudes::Vector{Float64} = SP.QUICK ? [0.002,0.005,0.01,0.02,0.05] : collect(0.002:0.002:0.05)
    signs::Vector{Int} = [-1, 1]
    delay_seeds::Vector{Int} = SP.QUICK ? [29] : [11, 17, 29]
    tail_cycles::Int = 6
    persistence_fraction::Float64 = 0.60
    baseline_margin::Float64 = 0.20
    bisection_steps::Int = SP.QUICK ? 4 : 7
    # Must be smaller than half of d_character-d_unstable (~0.028), so the
    # character neighbourhood cannot overlap the analytical separatrix.
    character_tolerance::Float64 = 0.02
    tail_std_max::Float64 = 0.12
    output_dir::String = joinpath(SP.FIGURES_DIR, "two_field_basin_boundary")
end

# Sub-sample peak times from a quadratic interpolation. This avoids quantising the
# measured basin coordinate at sample_dt or chorus-bin centres.
function refined_flashes(cfg, times, X, Y)
    signal = cfg.flash_variable == :x ? X : Y
    late = findall(times .>= cfg.t_transient)
    threshold = quantile(vec(signal[:,late]), cfg.flash_quantile)
    flashes = [Float64[] for _ in axes(signal,1)]
    for i in axes(signal,1), k in 2:size(signal,2)-1
        s0,s1,s2=signal[i,k-1],signal[i,k],signal[i,k+1]
        if times[k]>=cfg.t_transient && s1>threshold && s1>s0 && s1>=s2
            den=s0-2s1+s2
            δ=abs(den)>eps() ? clamp(0.5*(s0-s2)/den,-0.5,0.5) : 0.0
            push!(flashes[i],times[k]+δ*(times[k+1]-times[k]))
        end
    end
    flashes
end

function reference_chorus(cfg,times,flashes,reference)
    ref=[flashes[i] for i in reference]
    coarse,_,_=SP.collective_chorus(cfg,times,ref)
    refined=Float64[]
    for tc in coarse
        near=Float64[]
        for ft in ref
            j=searchsortedfirst(ft,tc)
            for q in (j-1,j)
                1<=q<=length(ft) && abs(ft[q]-tc)<=max(0.75,cfg.chorus_window) && push!(near,ft[q])
            end
        end
        length(near)>=max(3,round(Int,0.10length(reference))) && push!(refined,median(near))
    end
    refined
end

function delays_to_reference(flashes,chorus,ids)
    D=fill(NaN,length(ids),length(chorus))
    isempty(chorus) && return D
    halfwin=length(chorus)>2 ? 0.48median(diff(chorus)) : 2.0
    for (ii,i) in enumerate(ids), (n,tc) in enumerate(chorus)
        ft=flashes[i]; j=searchsortedfirst(ft,tc); best=NaN; ba=Inf
        for q in (j-1,j)
            if 1<=q<=length(ft) && abs(ft[q]-tc)<ba
                best=ft[q]-tc; ba=abs(best)
            end
        end
        ba<=halfwin && (D[ii,n]=best)
    end
    D
end

function analytical_geometry()
    path = joinpath(SP.FIGURES_DIR, "two_field_memory_validation", "homogeneous_equilibria.csv")
    isfile(path) || error("Run scripts/two_field_memory_validation.jl first; missing $path")
    a = readdlm(path, ',', Any, '\n'; header=true)[1]
    unstable = Float64[]; stablevals = Float64[]
    for row in eachrow(a)
        d = Float64(row[1]); isstable = lowercase(string(row[3])) == "true"
        if d > 0
            isstable ? push!(stablevals,d) : push!(unstable,d)
        end
    end
    isempty(unstable) && error("No positive unstable equilibrium found in $path")
    isempty(stablevals) && error("No positive stable character equilibrium found in $path")
    (;boundary=minimum(unstable),character=maximum(stablevals))
end

analytical_boundary() = analytical_geometry().boundary

function isotonic_probabilities(p,w)
    means=Float64[]; weights=Float64[]; starts=Int[]; stops=Int[]
    for i in eachindex(p)
        push!(means,p[i]); push!(weights,w[i]); push!(starts,i); push!(stops,i)
        while length(means)>=2 && means[end-1]>means[end]
            ww=weights[end-1]+weights[end]
            mm=(weights[end-1]*means[end-1]+weights[end]*means[end])/ww
            means[end-1]=mm; weights[end-1]=ww; stops[end-1]=stops[end]
            pop!(means); pop!(weights); pop!(starts); pop!(stops)
        end
    end
    out=zeros(Float64,length(p))
    for b in eachindex(means); out[starts[b]:stops[b]].=means[b]; end
    out
end

function kicked_simulation(cfg, case, lc, kick_amp, block_fraction, kick_sign, tail_cycles)
    N, dt = cfg.N, cfg.dt
    tend = cfg.t_end; nsteps = round(Int, tend/dt)+1
    every = max(1, round(Int, cfg.sample_dt/dt))
    phases = fill(0.0, N)
    x, y, _ = SP.init_phases_on_lc(cfg, lc, phases; seed=case.init_seed)
    coupling = SP.build_coupling_delay(cfg, case)
    yhist = zeros(N, nsteps); yhist[:,1] .= y
    ns = fld(nsteps-1, every)+1
    times=zeros(ns); X=zeros(N,ns); Y=zeros(N,ns); X[:,1].=x; Y[:,1].=y
    nk = clamp(round(Int, block_fraction*N), 1, N-1)
    block_first = fld(N-nk,2)+1; kicked = collect(block_first:block_first+nk-1)
    kickstep = round(Int, cfg.t_transient/dt)+1
    samp=1
    for step in 1:nsteps-1
        if step == kickstep
            for i in kicked
                phi = atan(y[i]-lc.cy, x[i]-lc.cx) + kick_sign*kick_amp
                x[i], y[i] = SP.sample_lc(lc, phi)
            end
            yhist[:,step] .= y
        end
        dx=similar(x); dy=similar(y)
        @inbounds for i in 1:N
            xi,yi=x[i],y[i]; coup=0.0
            for q in eachindex(coupling.neighbors[i])
                j=coupling.neighbors[i][q]
                yd=SP.delayed_value(yhist,step,coupling.taus[i][q],dt,j)
                coup += coupling.weights[i][q]*(yd-yi)
            end
            dx[i]=xi^2-cfg.beta*xi*yi
            dy[i]=cfg.eta*xi+cfg.rho*(xi-1)*yi+case.kappa*coup
        end
        x .+= dt.*dx; y .+= dt.*dy; yhist[:,step+1].=y
        if step % every == 0
            samp+=1; times[samp]=step*dt; X[:,samp].=x; Y[:,samp].=y
        end
    end
    times=times[1:samp]; X=X[:,1:samp]; Y=Y[:,1:samp]
    flashes=refined_flashes(cfg,times,X,Y)
    reference=setdiff(collect(1:N),kicked)
    chorus=reference_chorus(cfg,times,flashes,reference)
    D=delays_to_reference(flashes,chorus,kicked)
    post=findall(chorus .> cfg.t_transient + 0.25)
    first_delay=NaN; entry_delay=NaN; tail_delay=NaN; tail_std=Inf; tail_fraction=0.0
    if !isempty(post)
        v=filter(isfinite,D[:,first(post)])
        !isempty(v) && (first_delay=median(v))
        entrycols=length(post)>=2 ? post[2:min(length(post),4)] : Int[]
        if !isempty(entrycols)
            ev=filter(isfinite,vec(D[:,entrycols]))
            !isempty(ev) && (entry_delay=median(ev))
        end
        cols=post[max(1,end-tail_cycles+1):end]
        vals=D[:,cols]; good=filter(isfinite,vec(vals)); cyclemed=Float64[]
        for q in axes(vals,2)
            cv=filter(isfinite,vals[:,q]); !isempty(cv) && push!(cyclemed,median(cv))
        end
        if !isempty(good)
            tail_delay=median(good)
            tail_fraction=mean(abs.(good) .> max(cfg.chorus_window,SP.CHARACTER_TRACK_MIN_DELAY))
            tail_std=length(cyclemed)>1 ? std(cyclemed) : Inf
        end
    end
    post_cycle_delay=Float64[]; post_chorus_time=Float64[]
    for q in post
        cv=filter(isfinite,D[:,q])
        if !isempty(cv)
            push!(post_cycle_delay,median(cv)); push!(post_chorus_time,chorus[q])
        end
    end
    (; first_delay,entry_delay,tail_delay,tail_std,tail_fraction,nchorus=length(chorus),nk,kicked,
       post_cycle_delay,post_chorus_time)
end

function run(bc::BasinConfig=BasinConfig())
    mkpath(bc.output_dir)
    geometry=analytical_geometry(); boundary=geometry.boundary; character=geometry.character
    cfg=SP.SearchConfig(t_end=bc.settle_time+bc.observe_time,t_transient=bc.settle_time)
    lc=SP.build_limit_cycle(cfg)
    rows=NamedTuple[]
    controls=Dict{Tuple{Int,Float64},Any}()
    for seed in bc.delay_seeds, frac in bc.block_fractions
        case=(P=9,kappa=0.0485,tau0=1.185,sigtau=0.062,ell=12.0,delay_seed=seed,init_seed=202)
        controls[(seed,frac)]=kicked_simulation(cfg,case,lc,0.0,frac,1,bc.tail_cycles)
    end
    for seed in bc.delay_seeds, frac in bc.block_fractions, sgn in bc.signs, amp in bc.kick_amplitudes
        case=(P=9,kappa=0.0485,tau0=1.185,sigtau=0.062,ell=12.0,delay_seed=seed,init_seed=202)
        base=controls[(seed,frac)]
        r=kicked_simulation(cfg,case,lc,amp,frac,sgn,bc.tail_cycles)
        base_branch_distance=isfinite(base.tail_delay) ? min(abs(base.tail_delay-character),abs(base.tail_delay+character)) : Inf
        baseline_valid=isfinite(base.tail_delay) && abs(base.tail_delay)<boundary &&
                       base_branch_distance>bc.character_tolerance && base.tail_std<=bc.tail_std_max
        branch_distance=isfinite(r.tail_delay) ? min(abs(r.tail_delay-character),abs(r.tail_delay+character)) : Inf
        transitioned=baseline_valid && branch_distance<=bc.character_tolerance && r.tail_std<=bc.tail_std_max
        entry_perturbation=(isfinite(r.entry_delay)&&isfinite(base.entry_delay)) ? r.entry_delay-base.entry_delay : NaN
        push!(rows,(delay_seed=seed,block_fraction=frac,sign=sgn,kick_amplitude=amp,
                    first_delay=r.first_delay,entry_delay=r.entry_delay,entry_perturbation=entry_perturbation,
                    tail_delay=r.tail_delay,tail_std=r.tail_std,
                    branch_distance=branch_distance,tail_character_fraction=r.tail_fraction,
                    baseline_tail_delay=base.tail_delay,baseline_tail_character_fraction=base.tail_fraction,
                    baseline_valid=Int(baseline_valid),transitioned=Int(transitioned),chorus_events=r.nchorus,block_size=r.nk))
        println("seed=$seed width=$frac sign=$sgn kick=$amp first=$(r.first_delay) transitioned=$transitioned")
    end

    # Refine the first false->true transition in actuator space. The final bracket
    # is transformed to the physically comparable measured-delay coordinate.
    for seed in bc.delay_seeds, frac in bc.block_fractions, sgn in bc.signs
        base=controls[(seed,frac)]
        base_branch_distance=isfinite(base.tail_delay) ? min(abs(base.tail_delay-character),abs(base.tail_delay+character)) : Inf
        baseline_valid=isfinite(base.tail_delay) && abs(base.tail_delay)<boundary &&
                       base_branch_distance>bc.character_tolerance && base.tail_std<=bc.tail_std_max
        baseline_valid || continue
        g=sort(filter(r->r.delay_seed==seed&&r.block_fraction==frac&&r.sign==sgn,rows),by=r->r.kick_amplitude)
        it=findfirst(r->r.transitioned==1,g)
        isnothing(it) && continue
        hiamp=g[it].kick_amplitude
        loamp=it==1 ? 0.0 : g[it-1].kick_amplitude
        case=(P=9,kappa=0.0485,tau0=1.185,sigtau=0.062,ell=12.0,delay_seed=seed,init_seed=202)
        for _ in 1:bc.bisection_steps
            amp=(loamp+hiamp)/2
            r=kicked_simulation(cfg,case,lc,amp,frac,sgn,bc.tail_cycles)
            branch_distance=isfinite(r.tail_delay) ? min(abs(r.tail_delay-character),abs(r.tail_delay+character)) : Inf
            transitioned=branch_distance<=bc.character_tolerance && r.tail_std<=bc.tail_std_max
            entry_perturbation=(isfinite(r.entry_delay)&&isfinite(base.entry_delay)) ? r.entry_delay-base.entry_delay : NaN
            push!(rows,(delay_seed=seed,block_fraction=frac,sign=sgn,kick_amplitude=amp,
                        first_delay=r.first_delay,entry_delay=r.entry_delay,entry_perturbation=entry_perturbation,
                        tail_delay=r.tail_delay,tail_std=r.tail_std,
                        branch_distance=branch_distance,tail_character_fraction=r.tail_fraction,
                        baseline_tail_delay=base.tail_delay,baseline_tail_character_fraction=base.tail_fraction,
                        baseline_valid=1,transitioned=Int(transitioned),chorus_events=r.nchorus,block_size=r.nk))
            transitioned ? (hiamp=amp) : (loamp=amp)
            println("bisection seed=$seed width=$frac sign=$sgn kick=$amp measured=$(r.first_delay) transitioned=$transitioned")
        end
    end
    sort!(rows,by=r->(r.delay_seed,r.block_fraction,r.sign,r.kick_amplitude))
    names=propertynames(first(rows)); data=[getproperty(r,n) for r in rows,n in names]
    writedlm(joinpath(bc.output_dir,"basin_trials.csv"),vcat(permutedims(collect(names)),data),',')
    pdiag=scatter([abs(r.entry_perturbation) for r in rows],[r.tail_delay for r in rows],
              markercolor=[r.transitioned==1 ? :red : :steelblue for r in rows],
              group=[string("width=",r.block_fraction,", sign=",r.sign) for r in rows],
              xlabel="baseline-relative |entry perturbation|",ylabel="settled tail delay",
              title="Numerical basin test; analytical boundary = $(round(boundary,digits=3))",legend=:outerright)
    hline!(pdiag,[-character,character],color=:gray,linestyle=:dot,label="stable analytical characters")
    savefig(pdiag,joinpath(bc.output_dir,"entry_tail_diagnostic.png"))

    # Empirical reduced-coordinate calibration. Each spatial/sign family gets its
    # own character anchor, preventing waveform and block-size scaling from being
    # mistaken for a displacement of the reduced-model basin boundary.
    anchors=Dict{Tuple{Float64,Int},Float64}()
    for frac in bc.block_fractions, sgn in bc.signs
        v=[abs(r.entry_perturbation) for r in rows if r.block_fraction==frac && r.sign==sgn &&
           r.transitioned==1 && isfinite(r.entry_perturbation)]
        !isempty(v) && (anchors[(frac,sgn)]=median(v))
    end
    qrows=NamedTuple[]
    for r in rows
        anchor=get(anchors,(r.block_fraction,r.sign),NaN)
        q=isfinite(anchor)&&anchor>0&&isfinite(r.entry_perturbation) ? abs(r.entry_perturbation)/anchor : NaN
        push!(qrows,(delay_seed=r.delay_seed,block_fraction=r.block_fraction,sign=r.sign,
                    kick_amplitude=r.kick_amplitude,q=q,outcome=r.transitioned,anchor=anchor))
    end
    qn=propertynames(first(qrows)); qd=[getproperty(r,n) for r in qrows,n in qn]
    writedlm(joinpath(bc.output_dir,"committor_trials.csv"),vcat(permutedims(collect(qn)),qd),',')
    edges=collect(0.0:0.10:1.60); committor=NamedTuple[]
    for b in 1:length(edges)-1
        v=[r.outcome for r in qrows if isfinite(r.q)&&edges[b]<=r.q<edges[b+1]]
        !isempty(v) && push!(committor,(q_center=(edges[b]+edges[b+1])/2,count=length(v),probability=mean(v)))
    end
    # Weighted isotonic regression estimates a nondecreasing committor without
    # allowing a single noisy early success to force all later bins to one.
    mono=isotonic_probabilities([r.probability for r in committor],[r.count for r in committor])
    q50=NaN
    j=findfirst(>=(0.5),mono)
    !isnothing(j) && (q50=committor[j].q_center)
    cn=(:q_center,:count,:probability,:monotone_probability)
    cd=[Any[r.q_center r.count r.probability mono[i]] for (i,r) in enumerate(committor)]
    cdata=isempty(cd) ? Matrix{Any}(undef,0,4) : reduce(vcat,cd)
    writedlm(joinpath(bc.output_dir,"committor_bins.csv"),vcat(permutedims(collect(cn)),cdata),',')
    qu=boundary/character
    p=plot([r.q_center for r in committor],mono,marker=:circle,label="numerical committor",
           xlabel="normalized chorus-to-character progress q",ylabel="P(character)",
           title="Committor validation: q50=$(round(q50,digits=3)), analytical qu=$(round(qu,digits=3))",
           xlim=(0,1.6),ylim=(-0.03,1.03))
    scatter!(p,[r.q_center for r in committor],[r.probability for r in committor],label="raw bins",alpha=0.6)
    vline!(p,[qu],color=:black,linestyle=:dash,label="analytical qu=du/dc")
    hline!(p,[0.5],color=:gray,linestyle=:dot,label="50% committor")
    savefig(p,joinpath(bc.output_dir,"basin_boundary_validation.png"))
    thresholds=NamedTuple[]
    for seed in bc.delay_seeds, frac in bc.block_fractions, sgn in bc.signs
        g=sort(filter(r->r.delay_seed==seed&&r.block_fraction==frac&&r.sign==sgn,rows),by=r->r.kick_amplitude)
        first_transition=findfirst(r->r.transitioned==1,g)
        j=findfirst(q->g[q].transitioned==0&&g[q+1].transitioned==1,1:max(0,length(g)-1))
        if !isnothing(first_transition) && first_transition==1
            lo=0.0; hi=abs(g[1].entry_perturbation)
            estimate=(lo+hi)/2; uncertainty=abs(hi-lo)/2
        elseif isnothing(j)
            lo=hi=estimate=uncertainty=NaN
        else
            d0=abs(g[j].entry_perturbation); d1=abs(g[j+1].entry_perturbation)
            lo=min(d0,d1); hi=max(d0,d1)
            estimate=(lo+hi)/2; uncertainty=abs(hi-lo)/2
        end
        push!(thresholds,(delay_seed=seed,block_fraction=frac,sign=sgn,
              lower_delay=lo,upper_delay=hi,numerical_threshold=estimate,
              uncertainty=uncertainty,analytical_boundary=boundary))
    end
    tn=propertynames(first(thresholds)); td=[getproperty(r,n) for r in thresholds,n in tn]
    writedlm(joinpath(bc.output_dir,"basin_thresholds.csv"),vcat(permutedims(collect(tn)),td),',')
    open(joinpath(bc.output_dir,"committor_summary.txt"),"w") do io
        println(io,"analytical_q_boundary=",qu); println(io,"numerical_q50=",q50)
        println(io,"difference=",isfinite(q50) ? q50-qu : NaN)
    end
    (;outdir=bc.output_dir,boundary,character,analytical_q=qu,numerical_q50=q50,rows,thresholds)
end

end
