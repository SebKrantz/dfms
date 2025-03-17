#--------------------------------------------------------------------------
# PROCEDURES
#--------------------------------------------------------------------------

function InitCond(xNaN, r, p, optNaN, Rcon, q, NQ) # Rcon is R_mat, NQ is nQ

    rC = size(Rcon,2)
    pC = max(p,rC)

    x, indNaN = remNaNs_spline(xNaN, optNaN); 
        
    T, N = size(x);
    NM = N - NQ
    # Eigenval decomp of cov(x) = VDV', only r largest evals
    eig = eigen(cov(x));
    d = eig.values[1:r]; 
    v = eig.vectors[:, 1:r];
    
    # Observation equation   
    C = [v zeros(N,r*(pC-1))];
    
    # Static predictions
    f = x * v;

    ff = zeros(T-rC+1, 0);
    for kk = 0:rC-1
        ff = [ff f[rC-kk:end-kk,:]]; # Matrix of lagged factors
    end

    Rcon = kron(Rcon,eye(r));
    q = kron(q,ones(r,1));
    
    # This loops over the quarterly variables 
    for i = N-NQ+1:N
        xx_i = xNaN[rC:T, i]
        if sum(.!isnan.(xx_i)) < size(ff, 2)+2
            xx_i = x[rC:T, i]
        end
        ff_i = ff[.!isnan.(xx_i), :]
        xx_i = xx_i[.!isnan.(xx_i)] # Quarterly observations (no interpolation)
        iff_i = inv(ff_i'*ff_i)
        Cc = iff_i*ff_i'*xx_i # Coefficients from regressing quarterly observations on corresponding values of factors
        # This is restricted least squares with restrictions: Rcon * C_0 = q
        # The restrictions in Rcon (with -1 in the right places) make sense!
        Cc = Cc - iff_i*Rcon'*inv(Rcon*iff_i*Rcon')*(Rcon*Cc-q)
        C[i, 1:rC*r] = Cc' # This replaces the corresponding row. 
    end

    # This computes residuals based on the new C matrix (with and without missing values)
    res = x[rC:end, :] - ff * C[:, 1:rC*r]'
    resNaN = copy(res)
    resNaN[indNaN[rC:end, :]] .= NaN

    R = Diagonal(dropdims(nanvar(resNaN, dims = 1), dims = 1))
    R[NM+1:end, NM+1:end] .= 0
    C = [C [zeros(NM,rC*NQ); kron([1 2 3 2 1], eye(NQ))]]

    # Estimate A & Q from stacked F(t) = A*F(t-1) + e(t);
    z = f;
    Z = zeros(T-p, 0);
    for kk = 1:p
        Z = [Z z[p-kk+1:end-kk,:]]; # stacked regressors (lagged SPC)
    end
    z = z[p+1:end,:];
    ## run the var chi(t) = A*chi(t-1) + e(t);
    
    A = zeros(r*pC, r*pC);
    A_temp = inv(Z'*Z)*Z'*z;
    A[1:r,1:r*p] = A_temp';
    A[r+1:end,1:r*(pC-1)] = eye(r*(pC-1));

    temp = zeros(5*NQ, 5*NQ)
    temp[NQ+1:end,1:end-NQ] = eye(4*NQ)
    A = BlockDiagonal([A, temp])
    # Turn A into a normal matrix from BlockDiagonal
    A = Matrix(A)

    # Define Q matrix with dimensions pC*r+5*NQ x pC*r+5*NQ and fill it with zeros
    Q = zeros(pC*r+5*NQ, pC*r+5*NQ)
    # Calculate VAR residuals
    e = z - Z*A_temp
    # Extract covariance matrix of e and assign it to the top-left r x r block of Q
    Q[1:r, 1:r] = cov(e)
    # Extract diagonal elements of the nanvar of resNaN from column NM+1 to end and assign it to the bottom-right NQ x NQ block of Q
    Q[pC*r+1:pC*r+NQ, pC*r+1:pC*r+NQ] = Diagonal(dropdims(nanvar(resNaN[:, NM+1:end], dims = 1), dims = 1)) #  / 19

    # Calculate rp2
    # rp2 = (r*pC+5*NQ)^2

    # Initial conditions
    initZ = zeros(r*pC+5*NQ)
    initV = reshape(inv(I - kron(A, A)) * Q[:], r*pC+5*NQ, r*pC+5*NQ)

    return A, C, Q, R, initZ, initV
end

# y_est, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ
function EMstep(y, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ)
    n, T = size(y)
    rC = size(R_mat, 2)
    pC = max(p, rC)
    rp = r * p
    rpC = r * pC
    rC = rC * r

    # Running the Kalman filter with the current estimates of the parameters
    Zsmooth, Vsmooth, VVsmooth, loglik = runKF(y, A, C, Q, R, Z_0, V_0)
    
    EZZ   = Zsmooth[1:rpC,2:end] * Zsmooth[1:rpC,2:end]' + dropdims(sum(Vsmooth[1:rpC,1:rpC,2:end], dims = 3), dims = 3)
    EZZ_BB = Zsmooth[1:rpC,1:end-1] * Zsmooth[1:rpC,1:end-1]' + dropdims(sum(Vsmooth[1:rpC,1:rpC,1:end-1], dims = 3), dims = 3)
    EZZ_FB = Zsmooth[1:rpC,2:end] * Zsmooth[1:rpC,1:end-1]' + dropdims(sum(VVsmooth[1:rpC,1:rpC,:], dims = 3), dims = 3)
    
    A_new = copy(A)
    Q_new = copy(Q)
    
    A_new[1:r, 1:rp] = EZZ_FB[1:r, 1:rp] * inv(EZZ_BB[1:rp, 1:rp])
    Q_new[1:r, 1:r] = (EZZ[1:r, 1:r] - A_new[1:r, 1:rp] * EZZ_FB[1:r, 1:rp]') / T
    Q_new[rpC+1:rpC+nQ, rpC+1:rpC+nQ] = Diagonal(diag(Zsmooth[rpC+1:rpC+nQ, 2:end] * Zsmooth[rpC+1:rpC+nQ, 2:end]' + dropdims(sum(Vsmooth[rpC+1:rpC+nQ, rpC+1:rpC+nQ, 2:end], dims=3), dims = 3)) / T)
    
    Z_0 = Zsmooth[:, 1]
    V_0 = zeros(length(Z_0), length(Z_0))
    V_0[1:rpC, 1:rpC] = Vsmooth[1:rpC, 1:rpC, 1]
    V_0[rpC+1:end, rpC+1:end] = Diagonal(diag(Vsmooth[rpC+1:end, rpC+1:end, 1]))

    nanY = isnan.(y)
    y[nanY] .= 0

    C_new = C
    nM = n - nQ
    denom = zeros(nM*r, nM*r)
    nom = zeros(nM, r)

    for t=1:T
        nanYt = diagm(nanY[1:nM,t] .== 0)
        denom += kron(Zsmooth[1:r,t+1] * Zsmooth[1:r,t+1]' + Vsmooth[1:r,1:r,t+1], nanYt)
        nom += y[1:nM,t]*Zsmooth[1:r,t+1]';
    end
    
    vec_C = inv(denom) * nom[:]
    C_new = copy(C)
    C_new[1:nM, 1:r] = reshape(vec_C, nM, r)

    R_mat = kron(R_mat, eye(r))
    q = kron(q, ones(r, 1))
    for i=n-nQ+1:n
        denom = zeros(rC, rC)
        nom = zeros(1, rC)
    
        for t=1:T
            denom += !nanY[i,t] * (Zsmooth[1:rC, t+1]*Zsmooth[1:rC, t+1]' + Vsmooth[1:rC, 1:rC, t+1])
            nom += y[i, t]*Zsmooth[1:rC, t+1]'
        end
    
        C_i = inv(denom) * nom'
        C_i_constr = C_i - inv(denom)*R_mat'*inv(R_mat*inv(denom)*R_mat')*(R_mat*C_i-q)
        C_new[i, 1:rC] = C_i_constr
    end
    
    R_new = zeros(n, n)
    for t=1:T
        nanYt = diagm(nanY[:,t] .== 0)
        R_new += (y[:,t]-nanYt*C_new*Zsmooth[:,t+1]) * (y[:,t]-nanYt*C_new*Zsmooth[:,t+1])' +
                nanYt*C_new*Vsmooth[:,:,t+1]*C_new'*nanYt + (I(n)-nanYt)*R*(I(n)-nanYt)
    end
    
    R_new = Diagonal(diag(R_new) / T)
    R_new[n-nQ+1:end, n-nQ+1:end] .= 0

    return C_new, R_new, A_new, Q_new, loglik, Z_0, V_0
end

#--------------------------------------------------------------------------
# EM_CONVERGED Has EM converged?
# [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
#
# We have converged if the slope of the log-likelihood function falls below 'threshold', 
# i.e., |f(t) - f(t-1)| / avg < threshold,
# where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
# 'threshold' defaults to 1e-4.
#
# This stopping criterion is from Numerical Recipes in C p423
#
# If we are doing MAP estimation (using priors), the likelihood can decrase,
# even though the mode of the posterior is increasing.


function em_converged(loglik, previous_loglik, threshold=1e-4, check_increased=true)
    converged = false
    decrease = false

    if check_increased && (loglik - previous_loglik < -1e-3)
        # @printf("******likelihood decreased from %.4f to %.4f!\n", previous_loglik, loglik)
        decrease = true
    end
    
    delta_loglik = abs(loglik - previous_loglik)
    avg_loglik = (abs(loglik) + abs(previous_loglik) + eps())/2
    if (delta_loglik / avg_loglik) < threshold
        converged = true
    end
    
    return converged, decrease
end

     
        
