# -----------------------------
# Simulation
# -----------------------------
function simulate_case(cfg::SearchConfig, case, lc; keep_full=false, initial_phases=nothing)
    # lc is precomputed once in main() and reused for every case.
    # This avoids rebuilding the same single-firefly limit cycle during the parameter scan.
    x, y, init_phase = initial_phases === nothing ?
        init_random_on_lc(cfg, lc; seed=case.init_seed) :
        init_phases_on_lc(cfg, lc, initial_phases; seed=case.init_seed)
    coupling = build_coupling_delay(cfg, case)

    N = cfg.N
    dt = cfg.dt
    nsteps = Int(round(cfg.t_end/dt)) + 1
    sample_every = max(1, Int(round(cfg.sample_dt/dt)))
    nsamp = Int(floor((nsteps-1)/sample_every)) + 1

    # Full x-history is not needed for the delayed coupling; only y-history is stored.
    yhist = zeros(Float64, N, nsteps)
    yhist[:,1] .= y

    times = zeros(Float64, nsamp)
    Xs = zeros(Float64, N, nsamp)
    Ys = zeros(Float64, N, nsamp)
    samp = 1
    Xs[:,1] .= x; Ys[:,1] .= y

    for step in 1:(nsteps-1)
        # Explicit Euler method-of-steps. Use small dt for searches; refine candidates later.
        dx = similar(x); dy = similar(y)
        @inbounds for i in 1:N
            xi, yi = x[i], y[i]
            intrinsic_x = xi*xi - cfg.beta*xi*yi
            intrinsic_y = cfg.eta*xi + cfg.rho*(xi - 1.0)*yi
            coup = 0.0
            neigh = coupling.neighbors[i]
            wrow  = coupling.weights[i]
            taurow = coupling.taus[i]
            for q in eachindex(neigh)
                j = neigh[q]
                ydel = delayed_value(yhist, step, taurow[q], dt, j)
                coup += wrow[q] * (ydel - yi)
            end
            dx[i] = intrinsic_x
            dy[i] = intrinsic_y + case.kappa*coup
        end
        x .+= dt .* dx
        y .+= dt .* dy
        yhist[:,step+1] .= y
        if step % sample_every == 0
            samp += 1
            times[samp] = step*dt
            Xs[:,samp] .= x
            Ys[:,samp] .= y
        end
    end

    times = times[1:samp]
    Xs = Xs[:,1:samp]
    Ys = Ys[:,1:samp]
    return analyze_case(cfg, case, lc, times, Xs, Ys, coupling, init_phase; keep_full=keep_full)
end

# -----------------------------
