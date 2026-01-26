#--------------------------------------------------------------------------
# News_DFM_ML: News decomposition for DFM forecasts
# Implements Banbura & Modugno (2014) news decomposition
#--------------------------------------------------------------------------

using LinearAlgebra

# Only include if not already defined
if !@isdefined(para_const)
    include("para_const.jl")
end

"""
    News_DFM_ML(X_old, X_new, Q, t_fcst, v_news)

Compute news decomposition for forecast revision between two data vintages.

Arguments:
- X_old: T x N old data vintage (with more NaNs)
- X_new: T x N new data vintage (with releases)
- Q: Dict containing model parameters + metadata:
     C, R, Mx, Wx, Z_0, V_0, A, Q (state-space params)
     Groups: N-vector of group labels
     Series: N-vector of series names
- t_fcst: Target forecast time index
- v_news: Target variable index

Returns:
- y_old: Old forecast for target
- y_new: New forecast for target
- groupnews: News contribution by group
- singlenews: News contribution by variable
- gain: Gain weights for released variables
- gainSer: Series names for gain weights
- actual: Actual values of released variables
- fore: Old forecasts for released variables
- filt: Normalized gain weights
"""
function News_DFM_ML(X_old, X_new, Q, t_fcst, v_news)
    r = size(Q["C"], 2)
    T, N = size(X_new)
    gList = unique(Q["Groups"])
    groupnews = zeros(1, length(gList))
    singlenews = zeros(1, N)
    gain = Float64[]
    gainSer = String[]
    actual = zeros(N)
    fore = zeros(N)
    filt = zeros(N)

    if !isnan(X_new[t_fcst, v_news])
        # If new data contains target, news = observation - old forecast
        Res_old = para_const(X_old, Q, 0)
        temp = X_new[t_fcst, v_news] - Res_old["X_sm"][t_fcst, v_news]
        singlenews[v_news] = temp
        groupnews[findfirst(gList .== Q["Groups"][v_news])] = temp
        y_old = Res_old["X_sm"][t_fcst, v_news]
        y_new = X_new[t_fcst, v_news]
    else
        Mx = Q["Mx"]
        Wx = Q["Wx"]

        # Finding new releases
        miss_old = isnan.(X_old)
        miss_new = isnan.(X_new)
        temp = miss_old .- miss_new
        release_idx = findall(temp .== 1)
        t_miss = [idx[1] for idx in release_idx]
        v_miss = [idx[2] for idx in release_idx]

        lag = t_fcst .- t_miss
        k = isempty(lag) ? 0 : max(maximum(abs.(lag)), maximum(lag) - minimum(lag))

        C = Q["C"]
        R = Q["R"]'

        n_news = length(lag)

        # Run Kalman smoother on old and new data
        Res_old = para_const(X_old, Q, k)
        Res_new = para_const(X_new, Q, 0)

        y_old = Res_old["X_sm"][t_fcst, v_news]
        y_new = Res_new["X_sm"][t_fcst, v_news]

        if isempty(t_miss)
            return y_old, y_new, groupnews, singlenews, gain, gainSer, actual, fore, filt
        end

        # Computing innovations
        innov = zeros(n_news)
        for i = 1:n_news
            X_new_norm = (X_new[t_miss[i], v_miss[i]] - Mx[v_miss[i]]) / Wx[v_miss[i]]
            X_sm_norm = (Res_old["X_sm"][t_miss[i], v_miss[i]] - Mx[v_miss[i]]) / Wx[v_miss[i]]
            innov[i] = X_new_norm - X_sm_norm
        end

        # Computing weights for the news
        P = Res_old["P"][:, :, 2:end]

        # P1: covariance between target and release times
        P1 = zeros(r, n_news)
        for i = 1:n_news
            h = abs(t_fcst - t_miss[i])
            m = max(t_miss[i], t_fcst)
            if t_miss[i] > t_fcst
                Pp = P[1:r, h*r+1:h*r+r, m]'
            else
                Pp = P[1:r, h*r+1:h*r+r, m]
            end
            P1[:, i] = Pp * C[v_miss[i], 1:r]'
        end

        # P2: covariance matrix of releases
        P2 = zeros(n_news, n_news)
        WW = zeros(N, N)

        for i = 1:n_news
            for j = 1:n_news
                h = abs(lag[i] - lag[j])
                m = max(t_miss[i], t_miss[j])
                if t_miss[j] > t_miss[i]
                    Pp = P[1:r, h*r+1:(h+1)*r, m]'
                else
                    Pp = P[1:r, h*r+1:(h+1)*r, m]
                end
                if v_miss[i] == v_miss[j] && t_miss[i] != t_miss[j]
                    WW[v_miss[i], v_miss[j]] = 0
                else
                    WW[v_miss[i], v_miss[j]] = R[v_miss[i], v_miss[j]]
                end
                P2[i, j] = C[v_miss[i], 1:r]' * Pp * C[v_miss[j], 1:r] + WW[v_miss[i], v_miss[j]]
            end
        end

        # News and weights combined
        gain_vec = Wx[v_news] * C[v_news, 1:r]' * P1 * inv(P2)
        temp_news = gain_vec .* innov'

        # Aggregate news by variable and time
        singlenews_matrix = zeros(maximum(t_miss) - minimum(t_miss) + 1, N)
        actual = zeros(N)
        fore = zeros(N)
        filt = zeros(N)

        for i = 1:n_news
            row_idx = t_miss[i] - minimum(t_miss) + 1
            singlenews_matrix[row_idx, v_miss[i]] += temp_news[i]
            actual[v_miss[i]] = X_new[t_miss[i], v_miss[i]]
            fore[v_miss[i]] = Res_old["X_sm"][t_miss[i], v_miss[i]]
            filt[v_miss[i]] = gain_vec[i] / Wx[v_miss[i]]
        end

        singlenews = dropdims(sum(singlenews_matrix, dims=1), dims=1)

        # Aggregate by group
        for i = 1:length(gList)
            group_mask = Q["Groups"][v_miss] .== gList[i]
            if any(group_mask)
                groupnews[i] = sum(gain_vec[group_mask] .* innov[group_mask])
            end
        end

        # Get unique releases
        v_miss_unique, idx = unique_indices(v_miss)
        gain = gain_vec[idx]
        gainSer = Q["Series"][v_miss_unique]
    end

    return y_old, y_new, groupnews, singlenews, gain, gainSer, actual, fore, filt
end


"""
Helper function to get unique values and their first indices
"""
function unique_indices(v)
    seen = Dict{eltype(v), Int}()
    unique_vals = eltype(v)[]
    indices = Int[]

    for (i, val) in enumerate(v)
        if !haskey(seen, val)
            seen[val] = i
            push!(unique_vals, val)
            push!(indices, i)
        end
    end

    return unique_vals, indices
end
