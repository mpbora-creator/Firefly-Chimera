# -----------------------------
function make_positions(cfg::SearchConfig; seed=1)
    rng = MersenneTwister(seed)
    N = cfg.N
    if cfg.geometry == :ring
        angle = range(0, 2π, length=N+1)[1:N]
        return hcat(cos.(angle), sin.(angle))
    elseif cfg.geometry == :disk
        r = sqrt.(rand(rng, N))
        th = 2π .* rand(rng, N)
        return hcat(r .* cos.(th), r .* sin.(th))
    else
        error("Unknown geometry: $(cfg.geometry)")
    end
end

function build_coupling_delay(cfg::SearchConfig, case)
    # Optimized sparse neighbour-list coupling.
    # This avoids scanning all N oscillators for every i at every time step.
    # For ring_uniform, ell is deliberately ignored and should be singleton in the grid.
    N = cfg.N
    rng = MersenneTwister(case.delay_seed)

    neighbors = [Int[] for _ in 1:N]
    weights   = [Float64[] for _ in 1:N]
    taus      = [Float64[] for _ in 1:N]

    if cfg.kernel == :ring_uniform
        P = min(case.P, (N - 1) ÷ 2)
        for i in 1:N
            # neighbours on both sides of the ring
            for m in 1:P
                push!(neighbors[i], mod1(i - m, N))
                push!(neighbors[i], mod1(i + m, N))
            end
            w = 1.0 / length(neighbors[i])
            append!(weights[i], fill(w, length(neighbors[i])))
        end
    elseif cfg.kernel == :exponential
        cutoff = max(3 * case.P, 1)
        for i in 1:N
            raww = Float64[]
            for j in 1:N
                i == j && continue
                d = ringdist(i, j, N)
                if d <= cutoff
                    push!(neighbors[i], j)
                    push!(raww, exp(-d / case.ell))
                end
            end
            sw = sum(raww)
            append!(weights[i], sw > 0 ? raww ./ sw : raww)
        end
    else
        error("Unknown kernel: $(cfg.kernel)")
    end

    # Pairwise distributed delays. Delay heterogeneity is homogeneous: it does not define groups.
    for i in 1:N
        for j in neighbors[i]
            ξ = 2rand(rng) - 1
            dist_term = 0.015 * ringdist(i, j, N) / max(case.P, 1)
            push!(taus[i], max(cfg.dt, case.tau0 + dist_term + case.sigtau * ξ))
        end
    end

    return (; neighbors, weights, taus)
end

function delayed_value(hist::Matrix{Float64}, current_step::Int, delay::Float64, dt::Float64, i::Int)
    # Fractional-delay lookup with linear interpolation.
    # We need y_i(t - tau), where tau/dt is generally not an integer.
    # Let pos be the fractional history index corresponding to t - tau.
    # The returned value is linearly interpolated between the two nearest
    # stored history samples instead of rounding tau/dt to an integer step.
    pos = current_step - delay/dt

    # If the requested delayed time lies before the available history,
    # use the earliest stored state as a constant pre-history.
    if pos <= 1.0
        return hist[i, 1]
    end

    k0 = floor(Int, pos)
    k0 = clamp(k0, 1, current_step)
    k1 = clamp(k0 + 1, 1, current_step)
    α = pos - k0

    return (1.0 - α)*hist[i, k0] + α*hist[i, k1]
end
