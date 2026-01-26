#--------------------------------------------------------------------------
# KALMAN FILTER WITH LAG AUGMENTATION
# For computing cross-time covariances needed in news decomposition
#--------------------------------------------------------------------------

eye(n) = Matrix{Float64}(I, n, n)

"""
    runKF_lag(y, A, C, Q, R, x_0, Sig_0, k)

Kalman filter/smoother with state augmented to include k lags of factors.
This allows computing covariances P(t, t-h) for h = 0, ..., k.

Returns:
- xsmooth: Smoothed factor estimates (r x T+1)
- Vsmooth: Lag-augmented covariance matrices (r*(k+1) x r*(k+1) x T+1)
"""
function runKF_lag(y, A, C, Q, R, x_0, Sig_0, k)
    n, r_orig = size(C)

    if k > 0
        # Augment state to include k lags
        # New state: [f_t; f_{t-1}; ...; f_{t-k}]
        C = [C zeros(n, k * r_orig)]

        A_aug = zeros((k+1) * r_orig, (k+1) * r_orig)
        A_aug[1:r_orig, 1:r_orig] = A[1:r_orig, 1:min(r_orig, size(A, 2))]
        A_aug[r_orig+1:end, 1:k*r_orig] = eye(k * r_orig)
        A = A_aug

        Q_aug = zeros((k+1) * r_orig, (k+1) * r_orig)
        Q_aug[1:r_orig, 1:r_orig] = Q[1:r_orig, 1:r_orig]
        Q = Q_aug

        x_0 = [x_0[1:r_orig]; zeros(k * r_orig)]
        Sig_0_aug = zeros((k+1) * r_orig, (k+1) * r_orig)
        Sig_0_aug[1:r_orig, 1:r_orig] = Sig_0[1:r_orig, 1:r_orig]
        Sig_0 = Sig_0_aug
    end

    S = SKF_lag(y, C, R, A, Q, x_0, Sig_0)
    S = FIS_lag(y, C, R, A, Q, S)

    # Return only the first r components of smoothed state
    r = r_orig
    xsmooth = S.AmT[1:r, :]
    Vsmooth = S.PmT

    return xsmooth, Vsmooth
end


#______________________________________________________________________
# Standard Kalman Filter (forward pass)
#______________________________________________________________________
function SKF_lag(Y, Z, R, T, Q, A_0, P_0)
    n, m = size(Z)
    nobs = size(Y, 2)

    Am = fill(NaN, (m, nobs))
    Pm = fill(NaN, (m, m, nobs))
    AmU = fill(NaN, (m, nobs + 1))
    PmU = fill(NaN, (m, m, nobs + 1))
    loglik = 0.0

    Au = A_0
    Pu = P_0

    AmU[:, 1] = Au
    PmU[:, :, 1] = Pu

    y_t = Float64[]
    Z_t = zeros(0, m)
    PZF = zeros(m, 0)

    for t = 1:nobs
        # Prediction step
        A = T * Au
        P = T * Pu * T' + Q
        P = 0.5 * (P + P')

        # Handle missing data
        y_t, Z_t, R_t, L_t = MissData_lag(Y[:, t], Z, R)

        if isempty(y_t)
            Au = A
            Pu = P
        else
            PZ = P * Z_t'
            iF = pinv(Z_t * PZ + R_t)
            PZF = PZ * iF

            V = y_t - Z_t * A
            Au = A + PZF * V
            Pu = P - PZF * PZ'
            Pu = 0.5 * (Pu + Pu')

            # Use logabsdet for numerical stability (det can be slightly negative due to rounding)
            loglik += 0.5 * (logabsdet(iF)[1] - V' * iF * V)
        end

        Am[:, t] = A
        Pm[:, :, t] = P
        AmU[:, t+1] = Au
        PmU[:, :, t+1] = Pu
    end

    KZ = isempty(y_t) ? zeros(m, m) : PZF * Z_t

    return (Am = Am, Pm = Pm, AmU = AmU, PmU = PmU, loglik = loglik, KZ = KZ)
end


#______________________________________________________________________
# Fixed Interval Smoother (backward pass)
#______________________________________________________________________
function FIS_lag(Y, Z, R, T, Q, S)
    m, nobs = size(S.Am)

    AmT = zeros(m, nobs + 1)
    PmT = zeros(m, m, nobs + 1)
    PmT_1 = zeros(m, m, nobs)

    AmT[:, nobs+1] = S.AmU[:, nobs+1]
    PmT[:, :, nobs+1] = S.PmU[:, :, nobs+1]
    PmT_1[:, :, nobs] = (I(m) - S.KZ) * T * S.PmU[:, :, nobs]

    J_2 = S.PmU[:, :, nobs] * T' * pinv(S.Pm[:, :, nobs])

    for t = nobs:-1:1
        PmU = S.PmU[:, :, t]
        Pm1 = S.Pm[:, :, t]
        P_T = PmT[:, :, t+1]
        P_T1 = PmT_1[:, :, t]

        J_1 = J_2

        AmT[:, t] = S.AmU[:, t] + J_1 * (AmT[:, t+1] - T * S.AmU[:, t])
        PmT[:, :, t] = PmU + J_1 * (P_T - Pm1) * J_1'

        if t > 1
            J_2 = S.PmU[:, :, t-1] * T' * pinv(S.Pm[:, :, t-1])
            PmT_1[:, :, t-1] = PmU * J_2' + J_1 * (P_T1 - T * PmU) * J_2'
        end
    end

    return (AmT = AmT, PmT = PmT, PmT_1 = PmT_1, loglik = S.loglik)
end


#______________________________________________________________________
# Handle missing data
#______________________________________________________________________
function MissData_lag(y, C, R)
    ix = findall(.!isnan.(y))
    e = eye(length(y))
    L = e[:, ix]

    y = y[ix]
    C = C[ix, :]
    R = R[ix, ix]

    return y, C, R, L
end
