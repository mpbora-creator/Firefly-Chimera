# Argument-principle certification for the full N x N delayed network operator.

Base.@kwdef struct ArgumentPrincipleConfig
    real_left::Float64 = 1e-4
    real_right::Float64 = 1.0
    imag_max::Float64 = 20.0
    points_per_edge::Int = QUICK ? 60 : 220
    diffusion::Float64 = 0.0
    equilibrium_tolerance::Float64 = 1e-8
    boundary_logabs_floor::Float64 = -30.0
end

function delayed_weight_matrix(coupling,z)
    N=length(coupling.neighbors); W=zeros(ComplexF64,N,N)
    for i in 1:N,q in eachindex(coupling.neighbors[i])
        j=coupling.neighbors[i][q]; W[i,j]+=coupling.weights[i][q]*exp(-z*coupling.taus[i][q])
    end
    return W
end

function characteristic_matrix(z,state,p,coupling,L,D)
    q=-p.a .+3*p.b.*state.^2 .-5*p.c.*state.^4
    return z.*Matrix{ComplexF64}(I,length(state),length(state))-Diagonal(complex.(q))-D.*complex.(L)-p.Lambda.*delayed_weight_matrix(coupling,z)
end

function rectangle_contour(c)
    n=c.points_per_edge; a,b=c.real_left,c.real_right; w=c.imag_max
    return vcat(complex.(range(a,b,length=n),-w),complex.(b,range(-w,w,length=n)[2:end]),
                complex.(range(b,a,length=n)[2:end],w),complex.(a,range(w,-w,length=n)[2:end-1]))
end

function unwrap_angles(theta)
    out=copy(theta)
    for i in 2:length(out)
        delta=out[i]-out[i-1]
        while delta>pi; out[i]-=2pi; delta-=2pi; end
        while delta< -pi; out[i]+=2pi; delta+=2pi; end
    end
    return out
end

function argument_principle_count(state,p,coupling,L;acfg=ArgumentPrincipleConfig())
    contour=rectangle_contour(acfg); phases=Float64[]; logabs=Float64[]
    for z in contour
        la,sgn=logabsdet(characteristic_matrix(z,state,p,coupling,L,acfg.diffusion))
        push!(logabs,la); push!(phases,angle(sgn))
    end
    unwrapped=unwrap_angles(phases); winding=round(Int,(unwrapped[end]-unwrapped[1])/(2pi))
    return (;unstable_roots=winding,contour,phase=unwrapped,logabs,min_logabs=minimum(logabs))
end

function run_argument_principle_validation(;outdir=joinpath(FIGURES_DIR,"argument_principle_validation"),
        acfg=ArgumentPrincipleConfig())
    mkpath(outdir); h=run_heterogeneous_state_validation()
    if h.solved.residual>acfg.equilibrium_tolerance
        write_namedtuples_csv(joinpath(outdir,"root_certificate.csv"),[(;unstable_roots=missing,certified=false,reason="stationary_residual_too_large",
            real_left=acfg.real_left,real_right=acfg.real_right,imag_max=acfg.imag_max,min_logabs=NaN,
            stationary_residual=h.solved.residual,stationary_converged=false)])
        return (;heterogeneous=h,certificate=nothing,outdir)
    end
    cert=argument_principle_count(h.solved.state,h.p,h.result.coupling,h.solved.L;acfg)
    certified=cert.min_logabs>acfg.boundary_logabs_floor
    write_namedtuples_csv(joinpath(outdir,"root_certificate.csv"),[(;unstable_roots=certified ? cert.unstable_roots : missing,
        certified,reason=certified ? "ok" : "contour_too_close_to_root",
        real_left=acfg.real_left,real_right=acfg.real_right,imag_max=acfg.imag_max,min_logabs=cert.min_logabs,
        stationary_residual=h.solved.residual,stationary_converged=h.solved.converged)])
    rows=[(;real=real(cert.contour[i]),imag=imag(cert.contour[i]),unwrapped_phase=cert.phase[i],logabsdet=cert.logabs[i]) for i in eachindex(cert.contour)]
    write_namedtuples_csv(joinpath(outdir,"argument_contour.csv"),rows)
    p1=plot(real.(cert.contour),imag.(cert.contour);label=false,xlabel="Re λ",ylabel="Im λ",title="Certification contour")
    p2=plot(cert.phase;label=false,xlabel="contour index",ylabel="unwrapped arg det",title="Winding=$(cert.unstable_roots)")
    savefig(plot(p1,p2;layout=(1,2),size=(1000,420),margin=7Plots.mm),joinpath(outdir,"argument_principle_validation.png"))
    return (;heterogeneous=h,certificate=cert,certified,outdir)
end
