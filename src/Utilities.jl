# -----------------------------
# Utilities
# -----------------------------
wrap_angle(a) = atan(sin(a), cos(a))
ringdist(i, j, N) = min(abs(i-j), N-abs(i-j))

function robust_span(v; lo=0.05, hi=0.95)
    isempty(v) && return 0.0
    s = sort(collect(skipmissing(v)))
    isempty(s) && return 0.0
    n = length(s)
    ilo = clamp(Int(floor(lo*(n-1)))+1, 1, n)
    ihi = clamp(Int(ceil(hi*(n-1)))+1, 1, n)
    return s[ihi] - s[ilo]
end

function circular_std(theta)
    r = abs(mean(exp.(im .* theta)))
    return sqrt(max(0.0, -2log(max(r, 1e-12))))
end

function local_order(theta::AbstractVector, delta::Int)
    N = length(theta)
    Z = zeros(Float64, N)
    for i in 1:N
        s = 0.0 + 0.0im
        count = 0
        for m in -delta:delta
            j = mod1(i + m, N)
            s += exp(im*theta[j])
            count += 1
        end
        Z[i] = abs(s / count)
    end
    return Z
end

function phase_difference_profile(theta::AbstractVector)
    N = length(theta)
    return [wrap_angle(theta[mod1(i+1,N)] - theta[i]) for i in 1:N]
end

# -----------------------------
