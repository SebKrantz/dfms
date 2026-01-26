#--------------------------------------------------------------------------
# PROCEDURES for EM_DFM_SS_idio_restrMQ
# Mixed-frequency DFM with:
# - AR(1) idiosyncratic errors for monthly variables
# - MA(5) [1,2,3,2,1] aggregation for quarterly variables
#--------------------------------------------------------------------------

function InitCond(xNaN, r, p, optNaN, Rcon, q, nQ)

    rC = size(Rcon, 2)
    pC = max(p, rC)

    x, indNaN = remNaNs_spline(xNaN, optNaN)

    T, N = size(x)
    nM = N - nQ

    # Eigenval decomp of cov(x) = VDV', only r largest evals
    # Note: Julia eigen() returns eigenvalues in ascending order, so we need the last r
    eig = eigen(cov(x))
    d = eig.values[end-r+1:end]
    v = eig.vectors[:, end-r+1:end]

    # Observation equation
    C = [v zeros(N, r*(pC-1))]

    # Static predictions
    f = x * v

    ff = zeros(T-rC+1, 0)
    for kk = 0:rC-1
        ff = [ff f[rC-kk:end-kk, :]]
    end

    Rcon_kr = kron(Rcon, eye(r))
    q_kr = kron(q, ones(r, 1))

    # This loops over the quarterly variables
    for i = N-nQ+1:N
        xx_i = xNaN[rC:T, i]
        if sum(.!isnan.(xx_i)) < size(ff, 2) + 2
            xx_i = x[rC:T, i]
        end
        ff_i = ff[.!isnan.(xx_i), :]
        xx_i = xx_i[.!isnan.(xx_i)]
        iff_i = inv(ff_i' * ff_i)
        Cc = iff_i * ff_i' * xx_i
        # Restricted least squares
        Cc = Cc - iff_i * Rcon_kr' * inv(Rcon_kr * iff_i * Rcon_kr') * (Rcon_kr * Cc - q_kr)
        C[i, 1:rC*r] = Cc'
    end

    res = x[rC:end, :] - ff * C[:, 1:rC*r]'
    resNaN = copy(res)
    resNaN[indNaN[rC:end, :]] .= NaN

    R = Diagonal(dropdims(nanvar(resNaN, dims = 1), dims = 1))

    # Build C matrix: [C_factors | I_monthly | zeros | MA5_quarterly]
    eyeN = eye(N)
    eyeN = eyeN[:, 1:nM]  # Only monthly variables get identity
    C = [C eyeN]

    #--------------------------------------------------------------------------
    # Transition equation for monthly AR(1) idiosyncratic
    #--------------------------------------------------------------------------
    BM = zeros(nM, nM)
    SM = zeros(nM, nM)

    T_res, _ = size(resNaN)

    for i = 1:nM
        res_i = resNaN[:, i]
        # number of leading zeros
        leadZero_idx = findall(1:T_res .== cumsum(isnan.(res_i)))
        leadZero = isempty(leadZero_idx) ? 0 : maximum(leadZero_idx)
        endZero_idx = findall(1:T_res .== cumsum(isnan.(res_i[end:-1:1])))
        endZero = isempty(endZero_idx) ? 0 : maximum(endZero_idx)

        res_i_clean = res[max(1, leadZero+1):min(T_res, T_res-endZero), i]

        if length(res_i_clean) > 2
            BM[i, i] = (res_i_clean[1:end-1]' * res_i_clean[1:end-1]) \ (res_i_clean[1:end-1]' * res_i_clean[2:end])
            SM[i, i] = cov(res_i_clean[2:end] - res_i_clean[1:end-1] * BM[i, i])
        else
            BM[i, i] = 0.5
            SM[i, i] = 1.0
        end
    end

    initViM = Diagonal(1 ./ diag(eye(nM) - BM.^2)) .* SM

    # Add MA(5) structure for quarterly variables
    C = [C [zeros(nM, 5*nQ); kron(eye(nQ), [1 2 3 2 1])]]

    # Set observation noise to small value (errors in state)
    Rdiag = diag(R)
    sig_e = Rdiag[nM+1:N] / 19
    R = 1e-04 * eye(N)

    #--------------------------------------------------------------------------
    # Transition for quarterly MA(5) idiosyncratic
    #--------------------------------------------------------------------------
    rho0 = 0.1
    BQ = kron(eye(nQ), [[rho0 zeros(1, 4)]; [eye(4) zeros(4, 1)]])
    temp = zeros(5, 5)
    temp[1, 1] = 1
    SQ = kron(Diagonal((1 - rho0^2) * sig_e), temp)

    initViQ = reshape(inv(eye((5*nQ)^2) - kron(BQ, BQ)) * SQ[:], 5*nQ, 5*nQ)

    #--------------------------------------------------------------------------
    # Factor VAR: Estimate A & Q from stacked F(t) = A*F(t-1) + e(t)
    #--------------------------------------------------------------------------
    z = f
    Z = zeros(T - p, 0)
    for kk = 1:p
        Z = [Z z[p-kk+1:end-kk, :]]
    end
    z = z[p+1:end, :]

    A = zeros(r*pC, r*pC)
    A_temp = inv(Z' * Z) * Z' * z
    A[1:r, 1:r*p] = A_temp'
    if r*pC > r
        A[r+1:end, 1:r*(pC-1)] = eye(r*(pC-1))
    end

    Q = zeros(pC*r, pC*r)
    e = z - Z * A_temp
    Q[1:r, 1:r] = cov(e)

    # Initial conditions for factors
    rp2 = (r*pC)^2
    initV_fac = reshape(inv(eye(rp2) - kron(A, A)) * Q[:], r*pC, r*pC)

    # Block diagonal structure
    A = Matrix(BlockDiagonal([A, BM, BQ]))
    Q = Matrix(BlockDiagonal([Q, SM, SQ]))

    # Initial conditions
    initZ = zeros(size(A, 1))
    initV = Matrix(BlockDiagonal([initV_fac, initViM, initViQ]))

    return A, C, Q, R, initZ, initV
end


function EMstep(y, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ)

    n, T = size(y)
    rC = size(R_mat, 2)
    pC = max(p, rC)
    rp = r * p
    rpC = r * pC
    rC_r = rC * r
    nM = n - nQ

    # Running the Kalman filter with the current estimates of the parameters
    Zsmooth, Vsmooth, VVsmooth, loglik = runKF(y, A, C, Q, R, Z_0, V_0)

    # E-step: compute sufficient statistics for factor block
    EZZ = Zsmooth[1:rpC, 2:end] * Zsmooth[1:rpC, 2:end]' +
          dropdims(sum(Vsmooth[1:rpC, 1:rpC, 2:end], dims = 3), dims = 3)
    EZZ_BB = Zsmooth[1:rpC, 1:end-1] * Zsmooth[1:rpC, 1:end-1]' +
             dropdims(sum(Vsmooth[1:rpC, 1:rpC, 1:end-1], dims = 3), dims = 3)
    EZZ_FB = Zsmooth[1:rpC, 2:end] * Zsmooth[1:rpC, 1:end-1]' +
             dropdims(sum(VVsmooth[1:rpC, 1:rpC, :], dims = 3), dims = 3)

    # E-step: sufficient statistics for idiosyncratic block (diagonal)
    idio_start = rpC + 1
    idio_end = size(Zsmooth, 1)
    EZZ2 = diagm(diag(Zsmooth[idio_start:idio_end, 2:end] * Zsmooth[idio_start:idio_end, 2:end]')) +
           diagm(diag(dropdims(sum(Vsmooth[idio_start:idio_end, idio_start:idio_end, 2:end], dims = 3), dims = 3)))
    EZZ_BB2 = diagm(diag(Zsmooth[idio_start:idio_end, 1:end-1] * Zsmooth[idio_start:idio_end, 1:end-1]')) +
              diagm(diag(dropdims(sum(Vsmooth[idio_start:idio_end, idio_start:idio_end, 1:end-1], dims = 3), dims = 3)))
    EZZ_FB2 = diagm(diag(Zsmooth[idio_start:idio_end, 2:end] * Zsmooth[idio_start:idio_end, 1:end-1]')) +
              diagm(diag(dropdims(sum(VVsmooth[idio_start:idio_end, idio_start:idio_end, :], dims = 3), dims = 3)))

    # M-step: update A and Q for factors
    A_new = copy(A)
    Q_new = copy(Q)

    A_new[1:r, 1:rp] = EZZ_FB[1:r, 1:rp] * inv(EZZ_BB[1:rp, 1:rp])
    Q_new[1:r, 1:r] = (EZZ[1:r, 1:r] - A_new[1:r, 1:rp] * EZZ_FB[1:r, 1:rp]') / T

    # M-step: update A and Q for idiosyncratic (diagonal updates)
    A_new2 = EZZ_FB2 * diagm(1 ./ diag(EZZ_BB2))
    Q_new2 = (EZZ2 - A_new2 * EZZ_FB2') / T

    # Update monthly AR(1) parameters
    A_new[rpC+1:rpC+nM, rpC+1:rpC+nM] = A_new2[1:nM, 1:nM]
    Q_new[rpC+1:rpC+nM, rpC+1:rpC+nM] = Q_new2[1:nM, 1:nM]

    # Initial conditions
    Z_0_new = Zsmooth[:, 1]
    V_0_new = zeros(size(V_0))
    V_0_new[1:rpC, 1:rpC] = Vsmooth[1:rpC, 1:rpC, 1]
    V_0_new[rpC+1:end, rpC+1:end] = Diagonal(diag(Vsmooth[rpC+1:end, rpC+1:end, 1]))

    # Handle missing data
    nanY = isnan.(y)
    y_work = copy(y)
    y_work[nanY] .= 0

    C_new = copy(C)

    # Update C for monthly variables
    denom = zeros(nM*r, nM*r)
    nom = zeros(nM, r)

    for t = 1:T
        nanYt = diagm(Float64.(.!nanY[1:nM, t]))
        denom += kron(Zsmooth[1:r, t+1] * Zsmooth[1:r, t+1]' + Vsmooth[1:r, 1:r, t+1], nanYt)
        nom += y_work[1:nM, t] * Zsmooth[1:r, t+1]' -
               nanYt * (Zsmooth[rpC+1:rpC+nM, t+1] * Zsmooth[1:r, t+1]' + Vsmooth[rpC+1:rpC+nM, 1:r, t+1])
    end

    vec_C = inv(denom) * nom[:]
    C_new[1:nM, 1:r] = reshape(vec_C, nM, r)

    # Update C for quarterly variables with restrictions
    R_mat_kr = kron(R_mat, eye(r))
    q_kr = kron(q, ones(r, 1))

    for i = n-nQ+1:n
        denom_q = zeros(rC_r, rC_r)
        nom_q = zeros(1, rC_r)

        idx_jQ = i - nM
        i_idio_jQ = (rpC + nM + 5*(idx_jQ-1) + 1):(rpC + nM + 5*idx_jQ)

        # Update initial covariance for quarterly idiosyncratic
        V_0_new[i_idio_jQ, i_idio_jQ] = Vsmooth[i_idio_jQ, i_idio_jQ, 1]
        # Update AR and Q for quarterly idiosyncratic (first element of each 5-block)
        A_new[i_idio_jQ[1], i_idio_jQ[1]] = A_new2[i_idio_jQ[1] - rpC, i_idio_jQ[1] - rpC]
        Q_new[i_idio_jQ[1], i_idio_jQ[1]] = Q_new2[i_idio_jQ[1] - rpC, i_idio_jQ[1] - rpC]

        for t = 1:T
            nanYt_scalar = Float64(!nanY[i, t])
            denom_q += nanYt_scalar * (Zsmooth[1:rC_r, t+1] * Zsmooth[1:rC_r, t+1]' + Vsmooth[1:rC_r, 1:rC_r, t+1])
            nom_q += y_work[i, t] * Zsmooth[1:rC_r, t+1]'
            nom_q -= nanYt_scalar * ([1 2 3 2 1] * Zsmooth[i_idio_jQ, t+1] * Zsmooth[1:rC_r, t+1]' +
                                     [1 2 3 2 1] * Vsmooth[i_idio_jQ, 1:rC_r, t+1])
        end

        C_i = inv(denom_q) * nom_q'
        C_i_constr = C_i - inv(denom_q) * R_mat_kr' * inv(R_mat_kr * inv(denom_q) * R_mat_kr') * (R_mat_kr * C_i - q_kr)
        C_new[i, 1:rC_r] = vec(C_i_constr)
    end

    # R stays fixed (small observation noise)
    R_new = R

    return C_new, R_new, A_new, Q_new, Z_0_new, V_0_new, loglik
end


#--------------------------------------------------------------------------
# EM_CONVERGED Has EM converged?
#--------------------------------------------------------------------------
function em_converged(loglik, previous_loglik, threshold=1e-4, check_increased=true)
    converged = false
    decrease = false

    if check_increased && (loglik - previous_loglik < -1e-3)
        decrease = true
    end

    delta_loglik = abs(loglik - previous_loglik)
    avg_loglik = (abs(loglik) + abs(previous_loglik) + eps()) / 2
    if (delta_loglik / avg_loglik) < threshold
        converged = true
    end

    return converged, decrease
end
