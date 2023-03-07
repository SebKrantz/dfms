#--------------------------------------------------------------------------
# PROCEDURES
#--------------------------------------------------------------------------

function InitCond(xNaN, r, p, optNaN, i_idio) 

    x, indNaN = remNaNs_spline(xNaN, optNaN); 
        
    T, N = size(x);
    # Eigenval decomp of cov(x) = VDV', only r largest evals
    eig = eigen(cov(x));
    d = eig.values[1:r]; 
    v = eig.vectors[:, 1:r];
    
    # Observation equation   
    C = [v zeros(N,r*(p-1))];
    
    eyeN = eye(N)[:, i_idio .== 1]; # Removing columns as necessary
    C = [C eyeN];
    
    # Static predictions
    f = x * v;
    
    chi = x*v*v';
    res = x-chi;
    resNaN = res;
    resNaN[indNaN] .= NaN;
    
    R = Diagonal(dropdims(nanvar(resNaN, dims = 1), dims = 1));
    
    #--------------------------------------------------------------------------
    # Transition equation
    #--------------------------------------------------------------------------
    
    ii_idio = findall(i_idio .== 1);
    n_idio = length(ii_idio);
    B = zeros(n_idio, n_idio);
    S = zeros(n_idio, n_idio);
    
    for i = 1:n_idio;
        R[ii_idio[i], ii_idio[i]] = 1e-04;
    
        res_i = resNaN[:, ii_idio[i]];
        # number of leading zeros
        leadZero = maximum(findall(1:T .== cumsum(isnan.(res_i))), init = 0);
        endZero = maximum(findall(1:T .== cumsum(isnan.(res_i[end:-1:1]))), init = 0);
        
        res_i = res[leadZero+1:end-endZero, ii_idio[i]];
        
        B[i,i] = inv(res_i[1:end-1]'*res_i[1:end-1])*res_i[1:end-1]'*res_i[2:end];
        S[i,i] = cov(res_i[2:end]-res_i[1:end-1]*B[i,i]);
    end
     
    
    # Estimate A & Q from stacked F(t) = A*F(t-1) + e(t);
    z = f;
    Z = zeros(T-p, 0);
    for kk = 1:p
        Z = [Z z[p-kk+1:end-kk,:]]; # stacked regressors (lagged SPC)
    end
    z = z[p+1:end,:];
    ## run the var chi(t) = A*chi(t-1) + e(t);
    
    A = zeros(r*p, r*p);
    A_temp = inv(Z'*Z)*Z'*z;
    A[1:r,1:r*p] = A_temp';
    A[r+1:end,1:r*(p-1)] = eye(r*(p-1));
    
    Q = zeros(p*r, p*r);
    e = z  - Z*A_temp; # VAR residuals
    Q[1:r, 1:r] = cov(e); # VAR covariance matrix
    
    
    # Initial conditions
    initZ = zeros(size(Z,2)+n_idio, 1); # #[randn(1,r*(nlag+1))]';
    rp2 = (r*p)^2;
    initV1 = reshape(inv(eye(rp2) - kron(A,A))*Q[:], r*p, r*p);
    initV2 = Diagonal(1 ./ diag(eye(size(B,1)) - B.^2)) .* S;
    initV = BlockDiagonal([initV1, initV2]);
    A = BlockDiagonal([A, B]);
    Q = BlockDiagonal([Q, S]);

    return A, C, Q, R, initZ, initV
end



# y = y_est #, A, C, Q, R, Z_0, V_0, r, p, i_idio
function EMstep(y, A, C, Q, R, Z_0, V_0, r, p, i_idio)
    n, T = size(y)
    
    #Running the Kalman filter with the current estimates of the parameters
    Zsmooth, Vsmooth, VVsmooth, loglik = runKF(y, A, C, Q, R, Z_0, V_0)
    
    EZZ   = Zsmooth[1:r*p,2:end] * Zsmooth[1:r*p,2:end]' + dropdims(sum(Vsmooth[1:r*p,1:r*p,2:end], dims = 3), dims = 3)
    EZZ_BB = Zsmooth[1:r*p,1:end-1] * Zsmooth[1:r*p,1:end-1]' + dropdims(sum(Vsmooth[1:r*p,1:r*p,1:end-1], dims = 3), dims = 3)
    EZZ_FB = Zsmooth[1:r*p,2:end] * Zsmooth[1:r*p,1:end-1]' + dropdims(sum(VVsmooth[1:r*p,1:r*p,:], dims = 3), dims = 3)
    
    EZZ2    = diagm(diag(Zsmooth[r*p+1:end, 2:end] * Zsmooth[r*p+1:end, 2:end]')) +
                    diagm(diag(dropdims(sum(Vsmooth[r*p+1:end,r*p+1:end,2:end], dims = 3), dims = 3)))       
    EZZ_BB2 = diagm(diag(Zsmooth[r*p+1:end, 1:end-1] * Zsmooth[r*p+1:end, 1:end-1]')) +
                    diagm(diag(dropdims(sum(Vsmooth[r*p+1:end,r*p+1:end,1:end-1], dims = 3), dims = 3)))
    EZZ_FB2 = diagm(diag(Zsmooth[r*p+1:end, 2:end] * Zsmooth[r*p+1:end, 1:end-1]')) + 
                    diagm(diag(dropdims(sum(VVsmooth[r*p+1:end,r*p+1:end,:], dims = 3), dims = 3)))
    
    A_new = copy(A)
    Q_new = copy(Q)
    
    A_new[1:r, 1:r*p] = EZZ_FB[1:r,:] * inv(EZZ_BB)
    Q_new[1:r,1:r] = (EZZ[1:r,1:r] - A_new[1:r,1:r*p]*EZZ_FB[1:r,:]' ) / T

    A_new[r*p+1:end, r*p+1:end] = EZZ_FB2 * diagm(1 ./ diag(EZZ_BB2))
    Q_new[r*p+1:end, r*p+1:end] = (EZZ2 - A_new[r*p+1:end, r*p+1:end]*EZZ_FB2') / T

    
    nanY = isnan.(y)
    y[nanY] .= 0
    denom = zeros(n*r,n*r)
    nom = zeros(n,r)
    
    for t=1:T
        
        nanYt = diagm(nanY[:,t] .== 0)

        denom += kron(Zsmooth[1:r,t+1] * Zsmooth[1:r,t+1]' + Vsmooth[1:r,1:r,t+1], nanYt)

        nom += y[:,t] * Zsmooth[1:r,t+1]' -
            nanYt[:,i_idio] * (Zsmooth[r*p+1:end,t+1] * Zsmooth[1:r,t+1]' + Vsmooth[r*p+1:end,1:r,t+1])
    end
    
    vec_C = inv(denom) * nom[:]
    C_new = copy(C)
    C_new[1:n, 1:r] = reshape(vec_C, n, r)

    R_new = zeros(n, n)
    
    for t=1:T
        nanYt = diagm(nanY[:,t] .== 0)
        R_new += (y[:,t]-nanYt*C_new*Zsmooth[:,t+1]) * (y[:,t]-nanYt*C_new*Zsmooth[:,t+1])' +
                nanYt*C_new*Vsmooth[:,:,t+1]*C_new'*nanYt + (I(n)-nanYt)*R*(I(n)-nanYt)
    end
    
    R_new ./= T
    R_new = Diagonal(diag(R_new))

    # Initial conditions
    V_0 = Vsmooth[1:r*p,1:r*p,1]
    V2_0 = Diagonal(diag(Vsmooth[r*p+1:end,r*p+1:end,1]))
    V_0 = BlockDiagonal([V_0, V2_0])

    Z_0 = Zsmooth[:,1]

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

     
        
