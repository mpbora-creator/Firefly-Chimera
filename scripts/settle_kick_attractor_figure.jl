# Paired settle-kick trajectories for the calibrated two-field character levels.
ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"..","src","SarfatiPelegChimera.jl"))
include(joinpath(@__DIR__,"..","src","TwoFieldBasinBoundaryValidation.jl"))
using .TwoFieldBasinBoundaryValidation
using Plots, Statistics
const B=TwoFieldBasinBoundaryValidation
const S=Main.SarfatiPelegChimera

bc=B.BasinConfig(settle_time=260.0,observe_time=360.0,tail_cycles=6)
cfg=S.SearchConfig(t_end=bc.settle_time+bc.observe_time,t_transient=bc.settle_time)
lc=S.build_limit_cycle(cfg)
case=(P=9,kappa=0.0485,tau0=1.185,sigtau=0.062,ell=12.0,delay_seed=17,init_seed=202)
frac=0.45

control=B.kicked_simulation(cfg,case,lc,0.0,frac,1,bc.tail_cycles)
negative=B.kicked_simulation(cfg,case,lc,0.008,frac,-1,bc.tail_cycles)
positive=B.kicked_simulation(cfg,case,lc,0.00575,frac,1,bc.tail_cycles)
dc=B.analytical_geometry().character

cycle(v)=collect(1:length(v.post_cycle_delay))
p=plot(cycle(control),control.post_cycle_delay; marker=:circle,linewidth=2,
       color=:gray,label="zero-kick control",xlabel="chorus cycles after kick",
       ylabel="kicked-block median delay",title="Settle–kick convergence to delayed character states")
plot!(p,cycle(negative),negative.post_cycle_delay;marker=:circle,linewidth=2.5,
      color=:royalblue,label="negative kick Δφ = −0.008")
plot!(p,cycle(positive),positive.post_cycle_delay;marker=:circle,linewidth=2.5,
      color=:crimson,label="positive kick Δφ = +0.00575")
hline!(p,[dc];color=:crimson,linestyle=:dash,linewidth=2,label="analytical +d_c = +0.50167")
hline!(p,[-dc];color=:royalblue,linestyle=:dash,linewidth=2,label="analytical −d_c = −0.50167")
hline!(p,[0.0];color=:black,alpha=.35,linewidth=1,label=false)

outdir=joinpath(S.FIGURES_DIR,"settle_kick_attractor")
mkpath(outdir)
savefig(p,joinpath(outdir,"settle_kick_attractor_convergence.png"))
open(joinpath(outdir,"settle_kick_attractor_tracks.csv"),"w") do io
    println(io,"trial,cycle,chorus_time,block_median_delay")
    for (name,r) in (("control",control),("negative",negative),("positive",positive))
        for q in eachindex(r.post_cycle_delay)
            println(io,"$name,$q,$(r.post_chorus_time[q]),$(r.post_cycle_delay[q])")
        end
    end
end
println("control tail = ",control.tail_delay," ± ",control.tail_std)
println("negative tail = ",negative.tail_delay," ± ",negative.tail_std)
println("positive tail = ",positive.tail_delay," ± ",positive.tail_std)
println("analytical dc = ±",dc)
println("output = ",outdir)
