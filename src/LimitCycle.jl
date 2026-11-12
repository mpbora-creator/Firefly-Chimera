# -----------------------------
function f_single(x, y, beta, rho, eta)
    dx = x*x - beta*x*y
    dy = eta*x + rho*(x - 1.0)*y
    return dx, dy
end

function rk4_single_step(x, y, dt, beta, rho, eta)
    k1x,k1y = f_single(x,y,beta,rho,eta)
    k2x,k2y = f_single(x+0.5dt*k1x, y+0.5dt*k1y, beta,rho,eta)
    k3x,k3y = f_single(x+0.5dt*k2x, y+0.5dt*k2y, beta,rho,eta)
    k4x,k4y = f_single(x+dt*k3x, y+dt*k3y, beta,rho,eta)
    return x + dt*(k1x+2k2x+2k3x+k4x)/6, y + dt*(k1y+2k2y+2k3y+k4y)/6
end

function build_limit_cycle(cfg::SearchConfig)
    nsteps = Int(round(cfg.lc_relax_time/cfg.lc_dt))
    xs = Float64[]; ys = Float64[]
    x, y = 0.65, 1.0
    keep_start = max(1, nsteps - cfg.lc_keep_points + 1)
    for step in 1:nsteps
        x, y = rk4_single_step(x, y, cfg.lc_dt, cfg.beta, cfg.rho, cfg.eta)
        if step >= keep_start
            push!(xs, x); push!(ys, y)
        end
    end
    cx, cy = mean(xs), mean(ys)
    phi = atan.(ys .- cy, xs .- cx)
    order = sortperm(phi)
    return (; x=xs[order], y=ys[order], phi=phi[order], cx, cy)
end

function sample_lc(lc, phi)
    # Simple nearest-neighbour sampling on geometric phase. Robust enough for initialization.
    idx = searchsortedfirst(lc.phi, wrap_angle(phi))
    idx = clamp(idx, 1, length(lc.phi))
    return lc.x[idx], lc.y[idx]
end

function init_random_on_lc(cfg::SearchConfig, lc; seed=1)
    rng = MersenneTwister(seed)
    N = cfg.N
    x0 = zeros(Float64, N); y0 = zeros(Float64, N)
    phases = 2π .* rand(rng, N)
    for i in 1:N
        x0[i], y0[i] = sample_lc(lc, phases[i])
        x0[i] += 1e-4*randn(rng)
        y0[i] += 1e-4*randn(rng)
    end
    return x0, y0, phases
end

function init_phases_on_lc(cfg::SearchConfig, lc, phases::AbstractVector{<:Real}; seed=1)
    length(phases) == cfg.N || throw(ArgumentError("expected $(cfg.N) initial phases, got $(length(phases))"))
    rng = MersenneTwister(seed)
    x0 = zeros(Float64, cfg.N); y0 = zeros(Float64, cfg.N)
    phase0 = Float64.(phases)
    for i in 1:cfg.N
        x0[i], y0[i] = sample_lc(lc, phase0[i])
        x0[i] += 1e-4*randn(rng)
        y0[i] += 1e-4*randn(rng)
    end
    return x0, y0, phase0
end

# -----------------------------
