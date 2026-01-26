#--------------------------------------------------------------------------
# para_const: Run Kalman smoother with lag-augmented state
# Used for computing smoothed estimates and cross-time covariances
#--------------------------------------------------------------------------

using LinearAlgebra

# Only include if not already defined
if !@isdefined(runKF_lag)
    include("runKF_lag.jl")
end

"""
    para_const(X, P, lag)

Run Kalman filter/smoother on data X using parameters P, with state
augmented to include `lag` additional factor lags.

Arguments:
- X: T x N data matrix
- P: Named tuple or Dict containing model parameters:
     Z_0, V_0, A, C, Q, R, Mx, Wx
- lag: Number of additional lags to include in state (for cross-time covariances)

Returns:
- Res: Dict with:
  - P: Smoothed covariance matrices (with lags)
  - X_sm: Smoothed/fitted data in original scale
"""
function para_const(X, P, lag)
    # Extract parameters
    Z_0 = P["Z_0"]
    V_0 = P["V_0"]
    A = P["A"]
    C = P["C"]
    Q = P["Q"]
    R = P["R"]
    Mx = P["Mx"]
    Wx = P["Wx"]

    #--------------------------------------------------------------------------
    # Preparation of the data
    #--------------------------------------------------------------------------
    T, N = size(X)

    # Standardise x
    xNaN = (X .- Mx) ./ Wx

    y = xNaN'

    # Run Kalman filter with lag augmentation
    Zsmooth, P_smooth = runKF_lag(y, A, C, Q, R, Z_0, V_0, lag)

    Zsmooth = Zsmooth'
    x_sm = Zsmooth[2:end, :] * C'
    X_sm = Wx .* x_sm .+ Mx

    #--------------------------------------------------------------------------
    # Loading the results
    #--------------------------------------------------------------------------
    Res = Dict()
    Res["P"] = P_smooth
    Res["X_sm"] = X_sm

    return Res
end
