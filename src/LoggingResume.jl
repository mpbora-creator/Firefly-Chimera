# -----------------------------
# Output and plotting
# -----------------------------
function active_ell_values(cfg::SearchConfig, grid=active_grid())
    # For ring-uniform coupling, ell has no physical or numerical effect.
    # Returning a singleton removes duplicate cases from the grid.
    if cfg.kernel == :ring_uniform
        return [grid.ell_values[1]]
    elseif cfg.kernel == :exponential
        return grid.ell_values
    else
        error("Unknown kernel: $(cfg.kernel)")
    end
end

function build_cases(cfg::SearchConfig=CFG, grid=active_grid())
    cases = NamedTuple[]
    ell_list = active_ell_values(cfg, grid)
    for P in grid.P_values,
        kappa in grid.kappa_values,
        tau0 in grid.tau0_values,
        sigtau in grid.sigtau_values,
        ell in ell_list,
        delay_seed in grid.delay_seeds,
        init_seed in grid.init_seeds
        push!(cases, (; P, kappa, tau0, sigtau, ell, delay_seed, init_seed))
    end
    return cases
end


const CSV_HEADER = [:rank, :P, :kappa, :tau0, :sigtau, :ell, :delay_seed, :init_seed, :score,
    :target_chorus, :chorus_width, :synchrony_fraction,
    :chorus_score, :common_freq_score, :delay_persistence_score, :nonzero_delay_score,
    :di_track_score, :local_mixture_score, :char_activity_score,
    :two_cluster_penalty, :full_sync_penalty, :too_few_chorus_penalty,
    :chorus_count, :char_count, :concurrent_character_fraction,
    :stable_character_count, :stable_character_fraction,
    :stable_character_mean_delay, :stable_character_mean_std, :stable_character_mean_run,
    :mean_character_probability, :std_character_probability,
    :persistent_chorus_count, :switching_oscillator_count, :persistent_character_count,
    :mean_role_switch_count, :max_role_switch_count, :mean_role_lifetime,
    :delay_median, :delay_span, :R1_mean, :R2_mean, :Z_mean, :Z_std,
    :freq_mean, :freq_std, :freq_cv]

function summary_get(s, sym::Symbol, default)
    return hasproperty(s, sym) ? getproperty(s, sym) : default
end

function csv_values_from_summary(s; rank=0)
    return [rank, s.P, s.kappa, s.tau0, s.sigtau, s.ell, s.delay_seed, s.init_seed, s.score,
        s.target_chorus, s.chorus_width, s.synchrony_fraction,
        s.chorus_score, s.common_freq_score, s.delay_persistence_score, s.nonzero_delay_score,
        s.di_track_score, s.local_mixture_score, s.char_activity_score,
        s.two_cluster_penalty, s.full_sync_penalty, s.too_few_chorus_penalty,
        s.chorus_count, s.char_count, s.concurrent_character_fraction,
        s.stable_character_count, s.stable_character_fraction,
        s.stable_character_mean_delay, s.stable_character_mean_std, summary_get(s, :stable_character_mean_run, NaN),
        summary_get(s, :mean_character_probability, NaN), summary_get(s, :std_character_probability, NaN),
        summary_get(s, :persistent_chorus_count, 0), summary_get(s, :switching_oscillator_count, 0), summary_get(s, :persistent_character_count, 0),
        summary_get(s, :mean_role_switch_count, NaN), summary_get(s, :max_role_switch_count, 0), summary_get(s, :mean_role_lifetime, NaN),
        s.delay_median, s.delay_span, s.R1_mean, s.R2_mean, s.Z_mean, s.Z_std,
        s.freq_mean, s.freq_std, s.freq_cv]
end

function case_key(case)
    return string(case.P, "|", case.kappa, "|", case.tau0, "|", case.sigtau, "|",
                  case.ell, "|", case.delay_seed, "|", case.init_seed)
end

function case_key_from_summary(s)
    return string(s.P, "|", s.kappa, "|", s.tau0, "|", s.sigtau, "|",
                  s.ell, "|", s.delay_seed, "|", s.init_seed)
end

function ensure_live_csv(path=LIVE_CSV_PATH)
    mkpath(dirname(path))
    if !isfile(path) || filesize(path) == 0
        open(path, "w") do io
            println(io, join(string.(CSV_HEADER), ","))
            flush(io)
        end
    end
end

function append_live_csv(s, path=LIVE_CSV_PATH)
    ensure_live_csv(path)
    open(path, "a") do io
        println(io, join(csv_values_from_summary(s; rank=0), ","))
        flush(io)
    end
end

function parse_int_like(val::AbstractString)
    # Live CSV values may have been written as 10.0 instead of 10.
    # Accept both integer strings and integer-valued floating strings.
    txt = strip(replace(val, "\r" => ""))
    x = parse(Float64, txt)
    return Int(round(x))
end

function parse_csv_value(sym::Symbol, val::AbstractString)
    if sym in (:rank, :P, :delay_seed, :init_seed, :chorus_count, :char_count, :stable_character_count,
               :persistent_chorus_count, :switching_oscillator_count, :persistent_character_count, :max_role_switch_count)
        return parse_int_like(val)
    else
        return parse(Float64, strip(replace(val, "\r" => "")))
    end
end

function parse_key_value(sym::Symbol, val::AbstractString)
    txt = strip(replace(val, "\r" => ""))
    if sym in (:P, :delay_seed, :init_seed)
        return parse_int_like(txt)
    else
        return parse(Float64, txt)
    end
end

function load_logged_rows(path=LIVE_CSV_PATH)
    # Robust resume reader.  For deciding whether a case is complete, read only
    # the seven parameter-key columns.  Do not require all diagnostic columns to
    # parse, because a killed process may leave the last CSV row truncated.
    println("Resume live CSV path: ", abspath(path))
    if !isfile(path)
        println("Resume live CSV not found.")
        return Any[], Set{String}()
    end
    fsz = filesize(path)
    println("Resume live CSV size: ", fsz, " bytes")
    fsz == 0 && return Any[], Set{String}()

    lines = readlines(path)
    println("Resume live CSV lines including header: ", length(lines))
    length(lines) <= 1 && return Any[], Set{String}()

    header = Symbol.(strip.(split(strip(lines[1]), ",")))
    required = [:P, :kappa, :tau0, :sigtau, :ell, :delay_seed, :init_seed]
    col = Dict{Symbol,Int}()
    for (i,h) in enumerate(header)
        col[h] = i
    end
    missing_cols = [h for h in required if !haskey(col, h)]
    if !isempty(missing_cols)
        println("WARNING: live CSV header is missing resume columns: ", missing_cols)
        println("Header found: ", header)
        return Any[], Set{String}()
    end

    rows = Any[]
    completed = Set{String}()
    bad = 0
    duplicate = 0

    # Keep a lightweight summary when full columns are parseable; otherwise keep only the case.
    for (offset, line) in enumerate(lines[2:end])
        lineno = offset + 1
        txt = strip(line)
        isempty(txt) && continue
        parts = split(txt, ","; keepempty=true)
        if length(parts) < maximum(values(col))
            bad += 1
            continue
        end
        try
            case = (;
                P = parse_key_value(:P, parts[col[:P]]),
                kappa = parse_key_value(:kappa, parts[col[:kappa]]),
                tau0 = parse_key_value(:tau0, parts[col[:tau0]]),
                sigtau = parse_key_value(:sigtau, parts[col[:sigtau]]),
                ell = parse_key_value(:ell, parts[col[:ell]]),
                delay_seed = parse_key_value(:delay_seed, parts[col[:delay_seed]]),
                init_seed = parse_key_value(:init_seed, parts[col[:init_seed]]),
            )
            key = case_key(case)
            if key in completed
                duplicate += 1
            end
            push!(completed, key)

            # Try to reconstruct full summary for final summary/best-so-far.  If this fails,
            # still resume correctly because the completed key is already recorded.
            if length(parts) >= length(header)
                vals = Dict{Symbol,Any}()
                ok = true
                for (sym, val) in zip(header, parts)
                    try
                        vals[sym] = parse_csv_value(sym, val)
                    catch
                        ok = false
                        break
                    end
                end
                if ok
                    summary = NamedTuple{Tuple(header[2:end])}(Tuple(vals[h] for h in header[2:end]))
                    push!(rows, (; summary, case))
                end
            end
        catch err
            bad += 1
            # Print only a few examples to avoid flooding the terminal.
            if bad <= 3
                println("WARNING: could not read resume key on CSV line ", lineno, ": ", sprint(showerror, err))
                println("         line = ", txt[1:min(end, 160)])
            end
            continue
        end
    end

    println("Resume keys loaded: ", length(completed), " unique completed cases")
    if duplicate > 0
        println("Resume duplicate completed keys ignored: ", duplicate)
    end
    if bad > 0
        println("Resume malformed/truncated lines skipped: ", bad)
    end
    return rows, completed
end

function update_best_file(s, path=BEST_PATH)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "Best Sarfati–Peleg similarity case so far")
        println(io, @sprintf("score = %.8f", s.score))
        println(io, @sprintf("P = %d", s.P))
        println(io, @sprintf("kappa = %.8f", s.kappa))
        println(io, @sprintf("tau0 = %.8f", s.tau0))
        println(io, @sprintf("sigtau = %.8f", s.sigtau))
        println(io, @sprintf("ell = %.8f", s.ell))
        println(io, @sprintf("delay_seed = %d", s.delay_seed))
        println(io, @sprintf("init_seed = %d", s.init_seed))
        println(io, @sprintf("chorus_fraction = %.8f", s.synchrony_fraction))
        println(io, @sprintf("char_count = %d", s.char_count))
        println(io, @sprintf("freq_cv = %.8g", s.freq_cv))
        println(io, "")
        println(io, "Score components")
        println(io, @sprintf("chorus_score = %.8f", s.chorus_score))
        println(io, @sprintf("common_freq_score = %.8f", s.common_freq_score))
        println(io, @sprintf("delay_persistence_score = %.8f", s.delay_persistence_score))
        println(io, @sprintf("nonzero_delay_score = %.8f", s.nonzero_delay_score))
        println(io, @sprintf("di_track_score = %.8f", s.di_track_score))
        println(io, @sprintf("local_mixture_score = %.8f", s.local_mixture_score))
        println(io, @sprintf("char_activity_score = %.8f", s.char_activity_score))
        flush(io)
    end
end

function high_score_plot_path(s)
    mkpath(HIGH_SCORE_DIR)
    safe = @sprintf("score_%06.3f_P%d_k%.3f_tau%.3f_sig%.3f_dseed%d_iseed%d.png",
        s.score, s.P, s.kappa, s.tau0, s.sigtau, s.delay_seed, s.init_seed)
    return joinpath(HIGH_SCORE_DIR, safe)
end
