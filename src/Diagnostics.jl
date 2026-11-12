# ==========================================================
# Global-order plotting mode
# 0 = plot instantaneous R1(t), R2(t)
# 1 = plot time-averaged <R1>, <R2> as horizontal lines
# 2 = plot both instantaneous curves and their means
#
# This affects only the diagnostic figure display. It does not
# change theta, R1, R2, scoring, CSV output, or verification routines.
# ==========================================================
const GLOBAL_ORDER_PLOT_MODE = 1

# -----------------------------
function reconstruct_phase(Xs, Ys, lc)
    theta = atan.(Ys .- lc.cy, Xs .- lc.cx)
    return theta
end

function median_period(ft::AbstractVector{<:Real})
    if length(ft) >= 3
        ds = diff(collect(ft))
        good = ds[ds .> 0]
        return isempty(good) ? NaN : median(good)
    end
    return NaN
end

function detect_flashes(cfg::SearchConfig, times, Xs, Ys)
    N, M = size(Ys)
    signal = cfg.flash_variable == :x ? Xs : Ys
    late_cols = findall(times .>= cfg.t_transient)
    vals = vec(signal[:, late_cols])
    thr = quantile(vals, cfg.flash_quantile)
    flash_times = [Float64[] for _ in 1:N]
    flash_indices = [Int[] for _ in 1:N]

    for i in 1:N
        s = signal[i,:]
        for k in 2:(M-1)
            if times[k] >= cfg.t_transient && s[k] > thr && s[k] > s[k-1] && s[k] >= s[k+1]
                push!(flash_times[i], times[k])
                push!(flash_indices[i], k)
            end
        end
    end
    return flash_times, flash_indices, thr
end

function collective_chorus(cfg::SearchConfig, times, flash_times)
    allflashes = sort(reduce(vcat, flash_times))
    if isempty(allflashes)
        return Float64[], Float64[], Int[]
    end

    tmin, tmax = minimum(times), maximum(times)

    edges = collect(tmin:cfg.chorus_bin_dt:tmax)
    if edges[end] < tmax
        push!(edges, tmax)
    end

    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    counts = zeros(Int, length(centers))

    for t in allflashes
        if t < edges[1] || t > edges[end]
            continue
        end
        b = searchsortedlast(edges, t)
        if b == length(edges)
            b -= 1
        end
        if 1 <= b <= length(counts)
            counts[b] += 1
        end
    end

    minpeak = max(2, round(Int, cfg.min_chorus_peak_fraction * cfg.N))

    # Estimate natural period from individual oscillator flash intervals
    intervals = Float64[]
    for ft in flash_times
        if length(ft) >= 3
            append!(intervals, diff(ft))
        end
    end

    Tmed = isempty(intervals) ? 25.0 : median(intervals)

    # Only one chorus event allowed within about half a cycle
    min_sep_time = 0.45 * Tmed
    min_sep_bins = max(1, round(Int, min_sep_time / cfg.chorus_bin_dt))

    # Candidate local maxima
    cand = Int[]
    for k in 2:length(counts)-1
        if counts[k] >= minpeak &&
           counts[k] >= counts[k-1] &&
           counts[k] >= counts[k+1]
            push!(cand, k)
        end
    end

    if isempty(cand)
        return Float64[], centers, counts
    end

    # Non-maximum suppression:
    # choose strongest peak first, suppress nearby candidates
    order = sortperm(counts[cand], rev=true)
    selected = Int[]

    for q in order
        idx = cand[q]

        too_close = any(abs(idx - s) <= min_sep_bins for s in selected)

        if !too_close
            push!(selected, idx)
        end
    end

    sort!(selected)

    chorus = centers[selected]

    return chorus, centers, counts
end

function nearest_delay_to_chorus(t, chorus)
    isempty(chorus) && return missing
    idx = searchsortedfirst(chorus, t)
    candidates = Float64[]
    if idx <= length(chorus); push!(candidates, t - chorus[idx]); end
    if idx > 1; push!(candidates, t - chorus[idx-1]); end
    _, m = findmin(abs.(candidates))
    return candidates[m]
end

function classify_chorus_characters(cfg::SearchConfig, flash_times, chorus)
    chorus_delays = Float64[]
    character_delays = Float64[]
    character_counts_by_cycle = Dict{Int,Int}()
    chorus_count = 0
    char_count = 0
    for i in eachindex(flash_times)
        for t in flash_times[i]
            d = nearest_delay_to_chorus(t, chorus)
            ismissing(d) && continue
            if abs(d) <= cfg.chorus_window
                chorus_count += 1
                push!(chorus_delays, d)
            else
                char_count += 1
                push!(character_delays, d)
                k = searchsortedlast(chorus, t)
                character_counts_by_cycle[k] = get(character_counts_by_cycle, k, 0) + 1
            end
        end
    end
    concurrent_character_fraction = isempty(character_counts_by_cycle) ? 0.0 :
        mean(collect(values(character_counts_by_cycle)) .>= 2)
    total = chorus_count + char_count
    synchrony_fraction = total == 0 ? 0.0 : chorus_count/total
    return (; chorus_count, char_count, synchrony_fraction, chorus_delays,
            character_delays, concurrent_character_fraction)
end


function best_persistent_character_run(v::AbstractVector{<:Real}, cfg::SearchConfig)
    # Returns the best persistent off-chorus run for one oscillator.
    # v[n] is d_i(n), with NaN for missing cycles.
    # We allow a small number of missing/borderline cycles inside a run,
    # but the score is computed only from valid off-chorus delays.
    min_delay = max(cfg.chorus_window, CHARACTER_TRACK_MIN_DELAY)
    best_len = 0
    best_mean = 0.0
    best_std = Inf
    best_start = 0
    best_stop = 0

    run_vals = Float64[]
    run_start = 0
    last_valid = 0
    gap = 0

    function close_run!(stop_index)
        if length(run_vals) >= CHARACTER_TRACK_MIN_RUN
            μ = mean(run_vals)
            σ = length(run_vals) >= 2 ? std(run_vals) : 0.0
            # Require a nonzero phase delay and a reasonably horizontal d_i(n) band.
            if abs(μ) >= min_delay && σ <= CHARACTER_TRACK_STD_MAX
                # Prefer longer runs; use std only to break near ties.
                better = length(run_vals) > best_len ||
                         (length(run_vals) == best_len && σ < best_std)
                if better
                    best_len = length(run_vals)
                    best_mean = μ
                    best_std = σ
                    best_start = run_start
                    best_stop = stop_index
                end
            end
        end
    end

    for n in eachindex(v)
        d = Float64(v[n])
        is_character_cycle = isfinite(d) && abs(d) > min_delay
        if is_character_cycle
            if isempty(run_vals)
                run_start = n
            end
            push!(run_vals, d)
            last_valid = n
            gap = 0
        else
            if !isempty(run_vals)
                gap += 1
                if gap > CHARACTER_TRACK_MAX_GAP
                    close_run!(last_valid)
                    empty!(run_vals)
                    run_start = 0
                    last_valid = 0
                    gap = 0
                end
            end
        end
    end
    if !isempty(run_vals)
        close_run!(last_valid)
    end

    return (; length=best_len, mean_delay=best_mean, std_delay=best_std,
            start_cycle=best_start, stop_cycle=best_stop,
            is_stable=(best_len >= CHARACTER_TRACK_MIN_RUN && isfinite(best_std)))
end

function chorus_delay_tracks(cfg::SearchConfig, flash_times, chorus)
    # d_i(n): delay of oscillator i's nearest flash from chorus beat n.
    # Missing entries mean oscillator i did not flash close enough to that chorus cycle.
    N = length(flash_times)
    C = length(chorus)
    D = fill(NaN, N, C)
    if isempty(chorus)
        membership = character_membership_from_delay_matrix(cfg, D)
        return (; delay_matrix=D, cycles=collect(1:C), stable_delay_std=Float64[], stable_delay_mean=Float64[], stable_oscillators=Int[],
                stable_run_lengths=Int[], stable_start_cycles=Int[], stable_stop_cycles=Int[], membership...)
    end

    # Search half-window: half a chorus period, but not smaller than the chorus classification window.
    if C >= 3
        Tch = median(diff(chorus))
        halfwin = max(cfg.chorus_window, 0.48*Tch)
    else
        halfwin = max(cfg.chorus_window, 1.0)
    end

    for i in 1:N
        ft = flash_times[i]
        isempty(ft) && continue
        for (n, tc) in enumerate(chorus)
            idx = searchsortedfirst(ft, tc)
            best = NaN
            bestabs = Inf
            if idx <= length(ft)
                d = ft[idx] - tc
                if abs(d) < bestabs
                    best = d; bestabs = abs(d)
                end
            end
            if idx > 1
                d = ft[idx-1] - tc
                if abs(d) < bestabs
                    best = d; bestabs = abs(d)
                end
            end
            if bestabs <= halfwin
                D[i,n] = best
            end
        end
    end

    # Cycle-tracked stable delayed characters.
    # Instead of classifying each flash independently, follow every oscillator
    # across chorus cycles and look for horizontal nonzero d_i(n) bands.
    stable_std = Float64[]
    stable_mean = Float64[]
    stable_ids = Int[]
    stable_lengths = Int[]
    stable_start_cycles = Int[]
    stable_stop_cycles = Int[]
    for i in 1:N
        tr = best_persistent_character_run(view(D, i, :), cfg)
        if tr.is_stable
            push!(stable_ids, i)
            push!(stable_mean, tr.mean_delay)
            push!(stable_std, tr.std_delay)
            push!(stable_lengths, tr.length)
            push!(stable_start_cycles, tr.start_cycle)
            push!(stable_stop_cycles, tr.stop_cycle)
        end
    end

    membership = character_membership_from_delay_matrix(cfg, D)
    return (; delay_matrix=D, cycles=collect(1:C), stable_delay_std=stable_std,
            stable_delay_mean=stable_mean, stable_oscillators=stable_ids,
            stable_run_lengths=stable_lengths,
            stable_start_cycles=stable_start_cycles,
            stable_stop_cycles=stable_stop_cycles,
            membership...)
end


function character_membership_from_delay_matrix(cfg::SearchConfig, D::AbstractMatrix{<:Real})
    # Per-oscillator character-role diagnostic.
    # For every oscillator i, p_i is the fraction of observed chorus cycles
    # in which the oscillator is off-chorus with a nonzero delay.
    #
    # State convention for plotting:
    #   0.0  -> chorus-like cycle, |d_i(n)| <= chorus_window
    #   0.5  -> borderline/off-window but below character threshold
    #   1.0  -> character-like cycle, |d_i(n)| > max(chorus_window, CHARACTER_TRACK_MIN_DELAY)
    #   missing -> no flash assigned near that chorus cycle
    min_delay = max(cfg.chorus_window, CHARACTER_TRACK_MIN_DELAY)
    N, C = size(D)
    p = zeros(Float64, N)
    observed = zeros(Int, N)
    char_cycles = zeros(Int, N)
    chorus_cycles = zeros(Int, N)
    borderline_cycles = zeros(Int, N)
    state = Matrix{Union{Missing,Float64}}(missing, N, C)

    for i in 1:N, n in 1:C
        d = Float64(D[i,n])
        if !isfinite(d)
            state[i,n] = missing
            continue
        end
        observed[i] += 1
        if abs(d) > min_delay
            char_cycles[i] += 1
            state[i,n] = 1.0
        elseif abs(d) <= cfg.chorus_window
            chorus_cycles[i] += 1
            state[i,n] = 0.0
        else
            borderline_cycles[i] += 1
            state[i,n] = 0.5
        end
    end

    for i in 1:N
        p[i] = observed[i] == 0 ? NaN : char_cycles[i] / observed[i]
    end

    valid = observed .> 0
    persistent_chorus = [i for i in 1:N if valid[i] && p[i] <= 0.20]
    switchers = [i for i in 1:N if valid[i] && 0.20 < p[i] < 0.80]
    persistent_characters = [i for i in 1:N if valid[i] && p[i] >= 0.80]

    # ===== ROLE-SWITCH DYNAMICS DIAGNOSTIC =====
    # Count true role changes for each oscillator over chorus cycles.
    # We count transitions only between valid binary roles:
    #   0.0 = chorus-like, 1.0 = character-like.
    # Borderline cycles (0.5) and missing cycles are ignored when counting
    # switches. Thus p_i≈0.5 caused by a single long transition is separated
    # from p_i≈0.5 caused by rapid alternation.
    role_switch_count = zeros(Int, N)
    valid_role_cycle_count = zeros(Int, N)
    mean_role_lifetime = fill(NaN, N)
    max_role_lifetime = zeros(Int, N)

    for i in 1:N
        last_role = -1
        current_run = 0
        completed_runs = Int[]
        for n in 1:C
            st = state[i,n]
            role = -1
            if !(st isa Missing)
                if st == 0.0
                    role = 0
                elseif st == 1.0
                    role = 1
                end
            end
            role == -1 && continue
            valid_role_cycle_count[i] += 1
            if last_role == -1
                last_role = role
                current_run = 1
            elseif role == last_role
                current_run += 1
            else
                role_switch_count[i] += 1
                push!(completed_runs, current_run)
                last_role = role
                current_run = 1
            end
        end
        current_run > 0 && push!(completed_runs, current_run)
        if !isempty(completed_runs)
            mean_role_lifetime[i] = mean(completed_runs)
            max_role_lifetime[i] = maximum(completed_runs)
        end
    end
    valid_switch = valid_role_cycle_count .> 0
    mean_switch_count = any(valid_switch) ? mean(role_switch_count[valid_switch]) : 0.0
    max_switch_count = any(valid_switch) ? maximum(role_switch_count[valid_switch]) : 0
    mean_lifetime_all = any(isfinite.(mean_role_lifetime)) ? mean(mean_role_lifetime[isfinite.(mean_role_lifetime)]) : 0.0
    # ===== END ROLE-SWITCH DYNAMICS DIAGNOSTIC =====

    return (; character_probability=p,
            observed_cycle_count=observed,
            character_cycle_count=char_cycles,
            chorus_cycle_count=chorus_cycles,
            borderline_cycle_count=borderline_cycles,
            character_state_matrix=state,
            persistent_chorus=persistent_chorus,
            switchers=switchers,
            persistent_characters=persistent_characters,
            role_switch_count=role_switch_count,
            valid_role_cycle_count=valid_role_cycle_count,
            mean_role_lifetime=mean_role_lifetime,
            max_role_lifetime=max_role_lifetime,
            mean_switch_count=mean_switch_count,
            max_switch_count=max_switch_count,
            mean_role_lifetime_all=mean_lifetime_all)
end

function nan_to_missing_matrix(A::AbstractMatrix{<:Real})
    B = Matrix{Union{Missing,Float64}}(undef, size(A,1), size(A,2))
    for i in axes(A,1), j in axes(A,2)
        B[i,j] = isnan(A[i,j]) ? missing : Float64(A[i,j])
    end
    return B
end

function phase_frequency(theta, times, t0)
    late = findall(times .>= t0)
    length(late) < 4 && return zeros(size(theta,1))
    N = size(theta,1)
    ω = zeros(Float64, N)
    for i in 1:N
        ph = theta[i, late]
        un = zeros(length(ph)); un[1] = ph[1]
        for k in 2:length(ph)
            un[k] = un[k-1] + wrap_angle(ph[k] - ph[k-1])
        end
        ω[i] = (un[end]-un[1])/(times[late[end]]-times[late[1]])
    end
    return ω
end




# ============================================================
# Delay-persistence diagnostics
# ============================================================
# These diagnostics work directly with the cycle-resolved delay matrix
#
#     D[i,n] = d_i^(n) = t_i^(n) - t_chorus^(n),
#
# where i is oscillator index and n is chorus-cycle index.
#
# They are designed to quantify the Sarfati--Peleg role structure:
# persistent characters should have nearly horizontal delay tracks,
# switchers should have larger cycle-to-cycle delay variability, and
# groups with similar role dynamics should show correlated residuals.
#
# Returned quantities:
#   delay_variability[i]       = std_n D[i,n]
#   delay_residual_matrix[i,n] = D[i,n] - mean_n D[i,n]
#   delay_covariance_matrix[i,j] = cov_n(D[i,n],D[j,n])
#
function delay_persistence_diagnostics(D::AbstractMatrix{<:Real})
    N, C = size(D)

    delay_mean = fill(NaN, N)
    delay_variability = fill(NaN, N)
    delay_residual_matrix = fill(NaN, N, C)

    for i in 1:N
        vals = Float64[]
        for n in 1:C
            x = Float64(D[i,n])
            if isfinite(x)
                push!(vals, x)
            end
        end

        if !isempty(vals)
            μ = mean(vals)
            delay_mean[i] = μ
            delay_variability[i] = length(vals) >= 2 ? std(vals) : 0.0

            for n in 1:C
                x = Float64(D[i,n])
                if isfinite(x)
                    delay_residual_matrix[i,n] = x - μ
                end
            end
        end
    end

    delay_covariance_matrix = fill(NaN, N, N)
    for i in 1:N
        for j in i:N
            xs = Float64[]
            ys = Float64[]
            for n in 1:C
                xi = delay_residual_matrix[i,n]
                yj = delay_residual_matrix[j,n]
                if isfinite(xi) && isfinite(yj)
                    push!(xs, xi)
                    push!(ys, yj)
                end
            end

            if length(xs) >= 2
                cij = sum(xs .* ys) / (length(xs) - 1)
                delay_covariance_matrix[i,j] = cij
                delay_covariance_matrix[j,i] = cij
            elseif length(xs) == 1
                delay_covariance_matrix[i,j] = 0.0
                delay_covariance_matrix[j,i] = 0.0
            end
        end
    end

    return (; delay_mean, delay_variability, delay_residual_matrix, delay_covariance_matrix)
end



# Average delay covariance as a function of ring distance.
#
# For the covariance matrix C_ij, define the ring distance
#
#     r_ij = min(|i-j|, N-|i-j|).
#
# The returned profile is
#
#     C(r) = mean{ C_ij : r_ij = r }.
#
# r=0 is the diagonal and therefore represents the mean delay variance.
# r>0 measures how strongly delay fluctuations are correlated between
# oscillators separated by r sites on the ring.
function covariance_vs_ring_distance(C::AbstractMatrix{<:Real})
    N = size(C, 1)
    maxr = fld(N, 2)
    sums = zeros(Float64, maxr + 1)
    counts = zeros(Int, maxr + 1)

    for i in 1:N
        for j in 1:N
            cij = Float64(C[i,j])
            if !isfinite(cij)
                continue
            end
            raw = abs(i - j)
            r = min(raw, N - raw)
            if 0 <= r <= maxr
                sums[r+1] += cij
                counts[r+1] += 1
            end
        end
    end

    distances = collect(0:maxr)
    mean_cov = fill(NaN, maxr + 1)
    for k in eachindex(distances)
        if counts[k] > 0
            mean_cov[k] = sums[k] / counts[k]
        end
    end

    return (; distances, mean_cov, counts)
end

function analyze_case(cfg, case, lc, times, Xs, Ys, coupling, init_phase; keep_full=false)
    flash_times, flash_indices, flash_threshold = detect_flashes(cfg, times, Xs, Ys)

    # Use the original geometric/reconstructed phase for the chimera diagnostics.
    # This restores the earlier definition used in the parameter searches:
    #     theta_i(t) = atan2(y_i(t)-c_y, x_i(t)-c_x).
    # The flash/event phase routine is intentionally not used here, because it
    # changes the meaning of R1, R2, Z_i, phase snapshots, and the score.
    theta = reconstruct_phase(Xs, Ys, lc)

    chorus, chorus_bins, chorus_counts = collective_chorus(cfg, times, flash_times)
    classif = classify_chorus_characters(cfg, flash_times, chorus)
    delay_tracks = chorus_delay_tracks(cfg, flash_times, chorus)
    delay_diag = delay_persistence_diagnostics(delay_tracks.delay_matrix)
    delay_cov_distance = covariance_vs_ring_distance(delay_diag.delay_covariance_matrix)

    late_cols = findall(times .>= cfg.t_transient)
    late_theta = theta[:, late_cols]
    Mlate = length(late_cols)
    R1 = vec(abs.(mean(exp.(im .* late_theta), dims=1)))
    R2 = vec(abs.(mean(exp.(2im .* late_theta), dims=1)))
    Zmat = zeros(Float64, cfg.N, Mlate)
    Dmat = zeros(Float64, cfg.N, Mlate)
    for (kk,k) in enumerate(late_cols)
        Zmat[:,kk] .= local_order(theta[:,k], cfg.local_order_delta)
        Dmat[:,kk] .= phase_difference_profile(theta[:,k])
    end

    final_theta = theta[:, late_cols[end]]
    final_Z = Zmat[:, end]
    final_D = Dmat[:, end]
    ω = phase_frequency(theta, times, cfg.t_transient)
    freq_cv = std(ω) / max(abs(mean(ω)), 1e-9)

    # Sarfati–Peleg similarity score: not a proof, only a search guide.
    # Desired experimental analogue:
    #   1. a large-but-not-total chorus fraction near TARGET_CHORUS,
    #   2. common frequency across all oscillators,
    #   3. off-beat character flashes with persistent nonzero delays,
    #   4. cycle-resolved d_i(n) horizontal bands for delayed characters,
    #   5. spatial/local coherence variation without collapse to full synchrony.
    char_delays = classif.character_delays
    chorus_fraction = classif.synchrony_fraction
    total_flashes = classif.chorus_count + classif.char_count

    chorus_score = exp(-0.5*((chorus_fraction - TARGET_CHORUS)/CHORUS_WIDTH)^2)
    common_freq_score = exp(-0.5*(freq_cv/0.035)^2)

    # Character delays should be nonzero and not spread over the whole cycle.
    delay_span_val = isempty(char_delays) ? Inf : robust_span(char_delays)
    delay_median_val = isempty(char_delays) ? NaN : median(char_delays)
    delay_persistence_score = isempty(char_delays) ? 0.0 : exp(-0.5*(delay_span_val/0.85)^2)
    nonzero_delay_score = isempty(char_delays) ? 0.0 : clamp(abs(delay_median_val)/0.55, 0.0, 1.0)

    # d_i(n) score: reward stable off-beat oscillators whose delay is persistent over chorus cycles.
    stable_count = length(delay_tracks.stable_oscillators)
    stable_delay_std_mean = isempty(delay_tracks.stable_delay_std) ? Inf : mean(delay_tracks.stable_delay_std)
    stable_delay_abs_mean = isempty(delay_tracks.stable_delay_mean) ? 0.0 : mean(abs.(delay_tracks.stable_delay_mean))
    stable_run_mean = isempty(delay_tracks.stable_run_lengths) ? 0.0 : mean(delay_tracks.stable_run_lengths)
    role_p = delay_tracks.character_probability
    valid_role = isfinite.(role_p)
    mean_character_probability = any(valid_role) ? mean(role_p[valid_role]) : 0.0
    std_character_probability = count(valid_role) >= 2 ? std(role_p[valid_role]) : 0.0
    persistent_chorus_count = length(delay_tracks.persistent_chorus)
    switching_oscillator_count = length(delay_tracks.switchers)
    persistent_character_count = length(delay_tracks.persistent_characters)
    mean_role_switch_count = delay_tracks.mean_switch_count
    max_role_switch_count = delay_tracks.max_switch_count
    mean_role_lifetime = delay_tracks.mean_role_lifetime_all
    stable_fraction = stable_count / cfg.N
    stable_count_score = clamp(stable_fraction / 0.12, 0.0, 1.0)
    stable_std_score = isfinite(stable_delay_std_mean) ? exp(-0.5*(stable_delay_std_mean/0.18)^2) : 0.0
    stable_offset_score = clamp(stable_delay_abs_mean / 0.45, 0.0, 1.0)
    stable_run_score = clamp(stable_run_mean / 12.0, 0.0, 1.0)
    di_track_score = stable_count_score * stable_std_score * stable_offset_score * stable_run_score

    # Local mixture: prefer coexistence of highly coherent and lower-coherence regions.
    local_mixture_score = clamp(std(final_Z)/0.20, 0.0, 1.0) * exp(-0.5*((mean(final_Z)-0.82)/0.25)^2)

    # Ensure there are enough off-beat flashes to analyse, but do not reward excessive characters once present.
    char_activity_score = total_flashes == 0 ? 0.0 : clamp(classif.char_count / max(0.20*total_flashes, 1.0), 0.0, 1.0)

    # Penalties for states that are not Sarfati–Peleg-like.
    two_cluster_penalty = mean(R2) > 0.92 && mean(final_Z) > 0.92 ? 1.0 : 0.0
    full_sync_penalty = mean(R1) > 0.90 ? 1.0 : 0.0
    too_few_chorus_penalty = chorus_fraction < 0.20 ? 1.0 : 0.0

    score = W_CHORUS*chorus_score + W_FREQUENCY*common_freq_score +
            W_DELAY_SPAN*delay_persistence_score + W_DELAY_OFFSET*nonzero_delay_score +
            W_DI_TRACKS*di_track_score + W_LOCAL_MIX*local_mixture_score +
            W_CHAR_SYNC*classif.concurrent_character_fraction + W_ACTIVITY*char_activity_score -
            1.10*two_cluster_penalty - 1.10*full_sync_penalty - 0.80*too_few_chorus_penalty

    summary = (; case..., score, synchrony_fraction=chorus_fraction,
        target_chorus=TARGET_CHORUS, chorus_width=CHORUS_WIDTH,
        chorus_score, common_freq_score, delay_persistence_score, nonzero_delay_score,
        di_track_score, local_mixture_score, char_activity_score,
        two_cluster_penalty, full_sync_penalty, too_few_chorus_penalty,
        chorus_count=classif.chorus_count, char_count=classif.char_count,
        concurrent_character_fraction=classif.concurrent_character_fraction,
        stable_character_count=stable_count,
        stable_character_fraction=stable_fraction,
        stable_character_mean_delay=stable_delay_abs_mean,
        stable_character_mean_std=stable_delay_std_mean,
        stable_character_mean_run=stable_run_mean,
        mean_character_probability=mean_character_probability,
        std_character_probability=std_character_probability,
        persistent_chorus_count=persistent_chorus_count,
        switching_oscillator_count=switching_oscillator_count,
        persistent_character_count=persistent_character_count,
        mean_role_switch_count=mean_role_switch_count,
        max_role_switch_count=max_role_switch_count,
        mean_role_lifetime=mean_role_lifetime,
        delay_median=delay_median_val,
        delay_span=isfinite(delay_span_val) ? delay_span_val : NaN,
        R1_mean=mean(R1), R2_mean=mean(R2), Z_mean=mean(final_Z), Z_std=std(final_Z),
        delay_variability_mean=any(isfinite.(delay_diag.delay_variability)) ? mean(delay_diag.delay_variability[isfinite.(delay_diag.delay_variability)]) : NaN,
        delay_variability_std=any(isfinite.(delay_diag.delay_variability)) ? std(delay_diag.delay_variability[isfinite.(delay_diag.delay_variability)]) : NaN,
        delay_covariance_mean=any(isfinite.(delay_diag.delay_covariance_matrix)) ? mean(delay_diag.delay_covariance_matrix[isfinite.(delay_diag.delay_covariance_matrix)]) : NaN,
        delay_covariance_std=any(isfinite.(delay_diag.delay_covariance_matrix)) ? std(delay_diag.delay_covariance_matrix[isfinite.(delay_diag.delay_covariance_matrix)]) : NaN,
        delay_covariance_distance0=length(delay_cov_distance.mean_cov) >= 1 ? delay_cov_distance.mean_cov[1] : NaN,
        delay_covariance_distance1=length(delay_cov_distance.mean_cov) >= 2 ? delay_cov_distance.mean_cov[2] : NaN,
        delay_covariance_distanceP=(case.P <= length(delay_cov_distance.mean_cov)-1 ? delay_cov_distance.mean_cov[case.P+1] : NaN),
        freq_mean=mean(ω), freq_std=std(ω), freq_cv=freq_cv)

    if !keep_full
        return (; summary, case)
    end

    return (; summary, case, times, Xs, Ys, theta, flash_times, flash_indices, flash_threshold,
        chorus, chorus_bins, chorus_counts, classif, delay_tracks,
        delay_mean=delay_diag.delay_mean,
        delay_variability=delay_diag.delay_variability,
        delay_residual_matrix=delay_diag.delay_residual_matrix,
        delay_covariance_matrix=delay_diag.delay_covariance_matrix,
        delay_cov_distance=delay_cov_distance,
        R1, R2, Zmat, Dmat, final_theta, final_Z, final_D, omega=ω, coupling, lc, init_phase)
end

# ============================================================
# Global order-parameter verification diagnostics
# ============================================================
# These routines do NOT change any simulation, phase definition, score,
# or existing diagnostic.  They only analyse the already-computed geometric
# phase theta_i(t), R1(t), R2(t), and the collective flash-count signal.
#
# Purpose:
#   Test whether the periodic dips in R1,R2 are tied to the fast flash/reset
#   part of the relaxation cycle.  If the dips coincide with collective flash
#   peaks and with a transiently broad distribution of geometric phase, then
#   they should be interpreted as a geometric-phase/jump sensitivity effect,
#   not as loss of frequency locking.
# ============================================================

function nearest_series_value(xgrid::AbstractVector{<:Real}, y::AbstractVector, x::Real)
    isempty(xgrid) && return missing
    j = argmin(abs.(xgrid .- x))
    return y[j]
end

function circular_phase_histogram(ph::AbstractVector{<:Real}; nbins::Int=36)
    edges = collect(range(-pi, pi; length=nbins+1))
    counts = zeros(Int, nbins)
    for v in ph
        θ = wrap_angle(Float64(v))
        b = searchsortedlast(edges, θ)
        if b == nbins + 1
            b = nbins
        end
        if 1 <= b <= nbins
            counts[b] += 1
        end
    end
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    return (; centers, counts, edges)
end

function order_parameter_dip_diagnostics(res; nbins::Int=36)
    # Uses existing fields only:
    #   res.theta          geometric/reconstructed phase theta_i(t)
    #   res.R1,res.R2      order parameters already computed from theta
    #   res.chorus_bins    time grid for collective flash count
    #   res.chorus_counts  collective flash-count time series
    #   res.times          simulation sample times
    #   CFG.t_transient    start of analysed window
    late_cols = findall(res.times .>= CFG.t_transient)
    lt = res.times[late_cols]
    R1 = res.R1
    R2 = res.R2

    n = min(length(lt), length(R1), length(R2))
    lt = lt[1:n]
    R1 = R1[1:n]
    R2 = R2[1:n]

    flash_at_t = zeros(Float64, n)
    for k in 1:n
        v = nearest_series_value(res.chorus_bins, res.chorus_counts, lt[k])
        flash_at_t[k] = v isa Missing ? 0.0 : Float64(v)
    end

    V1 = 1.0 .- R1
    V2 = 1.0 .- R2

    corr_R1dip_flash = (std(V1) > 0 && std(flash_at_t) > 0) ? cor(V1, flash_at_t) : NaN
    corr_R2dip_flash = (std(V2) > 0 && std(flash_at_t) > 0) ? cor(V2, flash_at_t) : NaN

    # Representative snapshots: maximum R1, minimum R1, and recovery after minimum.
    imaxR = argmax(R1)
    iminR = argmin(R1)
    irecover = min(n, iminR + max(1, round(Int, 0.10*n)))

    snap_late_indices = unique([imaxR, iminR, irecover])
    snap_times = lt[snap_late_indices]
    snap_global_cols = late_cols[snap_late_indices]
    phase_hists = [circular_phase_histogram(res.theta[:, col]; nbins=nbins) for col in snap_global_cols]
    snap_labels = String[]
    for q in snap_late_indices
        if q == imaxR
            push!(snap_labels, "high R1")
        elseif q == iminR
            push!(snap_labels, "R1 dip")
        else
            push!(snap_labels, "after dip")
        end
    end

    # Chorus times and their corresponding R values.
    chorus_R1 = Float64[]
    chorus_R2 = Float64[]
    chorus_flash = Float64[]
    for tc in res.chorus
        j = argmin(abs.(lt .- tc))
        push!(chorus_R1, R1[j])
        push!(chorus_R2, R2[j])
        push!(chorus_flash, flash_at_t[j])
    end

    return (; times=lt, R1, R2, V1, V2, flash_count=flash_at_t,
            corr_R1dip_flash, corr_R2dip_flash,
            snap_times, snap_indices=snap_late_indices, snap_labels, phase_hists,
            chorus_times=res.chorus, chorus_R1, chorus_R2, chorus_flash)
end
