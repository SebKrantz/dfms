using LinearAlgebra

# Kalman Filter and Smoother

# X Data matrix (T x n)
# A Transition matrix (rp x rp)
# C Observation matrix (n x rp)
# Q State covariance (rp x rp)
# R Observation covariance (n x n)
# F_0 Initial state vector (rp x 1)
# P_0 Initial state covariance (rp x rp)
# retLL compute log likelihood

function SKFS(X::Matrix, A::Matrix, C::Matrix, Q::Matrix,
              R::Matrix, F_0::Vector, P_0::Matrix, retLL = false)

    T = size(X)[1]
    n = size(X)[2]
    rp = size(A)[1]
    n_c = Int64(0)

    # In internal code factors are Z (instead of F) and factor covariance V (instead of P),
    # to avoid confusion between the matrices and their predicted (p) and filtered (f) states.
    # Additionally the results matrices for all time periods have a T in the name.

    loglik = retLL ? 0 : NaN
    dn = 0
    if retLL
        dn = n * log(2.0*pi)
    end
    Zf = F_0
    Vf = P_0

    # Predicted state mean and covariance
    ZTp = zeros(rp,T)
    VTp = zeros(rp,rp,T)

    # Filtered state mean and covariance
    ZTf = zeros(rp,T)
    VTf = zeros(rp,rp,T)

    # Handling missing values in the filter
    arow = nmiss = findall(isfinite.(A[1,:]))
    if isempty(nmiss)
        error("Missing first row of transition matrix")
    end


    for i in 1:T

        # Run a prediction
        Zp = A * Zf
        Vp = A * Vf * transpose(A) + Q
        Vp += transpose(Vp) # Ensure symmetry
        Vp *= 0.5

        # If missing observations are present at some timepoints, exclude the
        # appropriate matrix slices from the filtering procedure.
        xt = X[i,:] 
        nmiss = findall(isfinite.(xt))
        n_c = length(nmiss)

        if n_c > 0
            if n_c == n
                Ci = C
                Ri = R
            else
                Ci = C[nmiss,arow]
                Ri = R[nmiss,nmiss]
                xt = xt[nmiss]
            end

            # Intermediate results
            VCt = Vp * transpose(Ci)
            S = inv(Symmetric(Ci * VCt + Ri))

            # Prediction error
            et = xt - Ci * Zp
            # Kalman gain
            K = VCt * S
            # Updated state estimate
            Zf = Zp + K * et
            # Updated state covariance estimate
            Vf = Vp - K * Ci * Vp
            Vf += transpose(Vf) # Ensure symmetry
            Vf *= 0.5

            # Compute likelihood. Skip this part if S is not positive definite.
            if retLL
                detS = det(S)
                if detS > 0
                    loglik += log(detS) - transpose(et) * S * et - dn # Standard Kalman Filter Likelihood
                end
            end

        else # If all missing: just prediction.
            Zf = Zp
            Vf = Vp
        end

        # Store predicted and filtered data needed for smoothing
        ZTp[:,i] = Zp
        VTp[:,:,i] = Vp
        ZTf[:,i] = Zf
        VTf[:,:,i] = Vf

    end

    if retLL
        loglik *= 0.5
    end

    # Smoothed state mean and covariance
    ZsT = zeros(rp,T)
    VsT = zeros(rp,rp,T)
    VVsT = zeros(rp,rp,T) # Cov(Z_t, Z_t-1), used in EM

    # Initialize smoothed data with last observation of filtered data
    ZsT[:,T] = Zf
    VsT[:,:,T] = Vf
    K = (n_c == 0) ? zeros(rp,rp) : K * Ci
    VVsT[:,:,T] = (I(rp) - K) * A * VTf[:,:,T-1]

    At = transpose(A)
    Jimt = transpose(VTf[:,:,T-1] * At * inv(VTp[:,:,T]))

    # Smoothed state variable and covariance
    for i = T-1:-1:1
        Vf = VTf[:,:,i]
        Ji = transpose(Jimt)
        ZsT[:,i] = ZTf[:,i] + Ji * (ZsT[:,i+1] - ZTp[:,i+1])
        VsT[:,:,i] = Vf + Ji * (VsT[:,:,i+1] - VTp[:,:,i+1]) * Jimt
        # Cov(Z_t, Z_t-1): Needed for EM
        if(i > 1)
            Jimt = transpose(VTf[:,:,i-1] * At * inv(VTp[:,:,i]))
            VVsT[:,:,i] = Vf * Jimt + Ji * (VVsT[:,:,i+1] - A * Vf) * Jimt
        end
    end

    # Smoothed value for t = 0
    Vp = VTp[:,:,1]
    Jimt = transpose(P_0 * At * inv(Vp))
    VVsT[:,:,1] = Vf * Jimt + Ji * (VVsT[:,:,2] - A * Vf) * Jimt
    Ji = transpose(Jimt)

    return (F = transpose(ZTf),
            P = VTf,
            F_pred = transpose(ZTp),
            P_pred = VTp,
            F_smooth = transpose(ZsT),
            P_smooth = VsT,
            PPm_smooth = VVsT,
            F_smooth_0 = F_0 + Ji * (ZsT[:,1] - ZTp[:,1]),
            P_smooth_0 = P_0 + Ji * (VsT[:,:,1] - Vp) * Jimt,
            loglik = loglik)
end
