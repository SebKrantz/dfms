#--------------------------------------------------------------------------
# KALMAN FILTER
#--------------------------------------------------------------------------
# runKF(y, A, C, Q, R, Z_0, V_0)
function runKF(y, A, C, Q, R, Z_0, V_0)

    S = SKF(y, C, R, A, Q, Z_0, V_0)
    S = FIS(y, C, R, A, Q, S)

    return (xsmooth = S.ZmT,
            Vsmooth = S.PmT,
            VVsmooth = S.PmT_1,
            loglik = S.loglik)
end 


#______________________________________________________________________
# Kalman filter for stationary systems with time-varying system matrices
# and missing data.
#
# The model is        y_t   = C * z_t + eps_t       
#                     z_t+1 = A * z_t + u_t       
#
#______________________________________________________________________
# INPUT  
#        Y         Data                                 (nobs x n)  
# OUTPUT 
#        S.Zm       Predicted state vector  Z_t|t-1      (nobs x m)  
#        S.ZmU      Filtered  state vector  Z_t|t        (nobs+1 x m)  
#        S.Pm       Predicted covariance of Z_t|t-1      (nobs x m x m)  
#        S.PmU      Filtered  covariance of Z_t|t        (nobs+1 x m x m)  
#        S.loglik   Value of likelihood function

# SKF(y, C, R, A, Q, Z_0, V_0)
function SKF(Y, C, R, A, Q, Z_0, P_0)
    n, m = size(C)
    nobs  = size(Y, 2)

    Zm = fill(NaN, (m, nobs))
    Pm = fill(NaN, (m, m, nobs))
    ZmU = fill(NaN, (m, nobs+1))
    PmU = fill(NaN, (m, m, nobs+1))
    loglik = 0

    Zu = Z_0
    Pu = P_0

    ZmU[:,1]   .= Zu
    PmU[:,:,1] .= Pu

    y_t = 0; C_t = 0; PCF = 0

    for t in 1:nobs
        Z = A * Zu
        P = A * Pu * transpose(A) + Q
        P = (P + transpose(P)) / 2
        
        y_t, C_t, R_t, L_t  = MissData(Y[:,t], C, R)

        if isempty(y_t)
            Zu = Z
            Pu = P
        else
            PC  = P * transpose(C_t)
            iF  = inv(C_t * PC + R_t)
            PCF = PC * iF

            V   = y_t - C_t * Z
            Zu  = Z  + PCF * V
            Pu  = P  - PCF * PC'
            Pu  = (Pu + transpose(Pu)) / 2
            loglik += 0.5 * (log(det(iF)) - first(transpose(V) * iF * V))
        end

        Zm[:,t]   .= Z
        Pm[:,:,t] .= P
        ZmU[:,t+1]   .= Zu
        PmU[:,:,t+1] .= Pu
    end 

    KC = isempty(y_t) ? zeros(m,m) : PCF * C_t
    return (Zm=Zm, Pm=Pm, ZmU=ZmU, PmU=PmU, loglik=loglik, KC=KC)
end

 

#______________________________________________________________________
# Fixed intervall smoother (see Harvey, 1989, p. 154)
# FIS returns the smoothed state vector ZmT and its covar matrix PmT             
# Use this in conjnuction with function SKF
#______________________________________________________________________
# INPUT  
#        Y         Data                                 (nobs x n)  
#        S Estimates from Kalman filter SKF                                                          
#          S.Zm   : Estimates     z_t|t-1                  (nobs x m) 
#          S.Pm   : P_t|t-1 = Cov(z_t|t-1)             (nobs x m x m)
#          S.ZmU  : Estimates     z_t|t                    (nobs x m) 
#          S.PmU  : P_t|t   = Cov(z_t|t)               (nobs x m x m)       
# OUTPUT 
#        S Smoothed estimates added to above
#          S.ZmT  : Estimates     z_t|T                    (nobs x m) 
#          S.PmT :  P_t|T   = Cov(z_t|T)               (nobs x m x m)
#          S.PmT_1 : Cov(z_tz_t-1|T)
#        where m is the dim of state vector and t = 1 ...T is time

# FIS(y, C, R, A, Q, S)
function FIS(Y, C, R, A, Q, S)
    m, nobs = size(S.Zm)
    ZmT = zeros(m, nobs+1)
    PmT = zeros(m, m, nobs+1)
    PmT_1 = zeros(m, m, nobs)
    ZmT[:, nobs+1] = vec(S.ZmU[:, nobs+1])
    PmT[:, :, nobs+1] = copy(S.PmU[:, :, nobs+1])
    PmT_1[:, :, nobs] = (I(m) - S.KC) * A * S.PmU[:, :, nobs]

    J_2 = S.PmU[:, :, nobs] * A' * pinv(S.Pm[:, :, nobs])

    for t = nobs:-1:1
        PmU = S.PmU[:, :, t]
        Pm1 = S.Pm[:, :, t]
        P_T = PmT[:, :, t+1]
        P_T1 = PmT_1[:, :, t]

        J_1 = J_2

        ZmT[:, t] = S.ZmU[:, t] + J_1 * (ZmT[:, t+1] - A * S.ZmU[:, t])
        PmT[:, :, t] = PmU + J_1 * (P_T - Pm1) * J_1'

        if t > 1
            J_2 = S.PmU[:, :, t-1] * A' * pinv(S.Pm[:, :, t-1])
            PmT_1[:, :, t-1] = PmU * J_2' + J_1 * (P_T1 - A * PmU) * J_2'
        end
    end

    return (ZmT = ZmT, PmT = PmT, PmT_1 = PmT_1, loglik = S.loglik)
end
# The code is very similar to the Matlab version, with a few minor changes:

#     The squeeze function is not needed in Julia since it supports slicing of arbitrary dimensions.
#     The vec function is used to flatten a 2D array into a 1D vector.
#     The copy function is used to make a deep copy of an array, since the = operator in Julia creates a shallow copy.
#     The eye function is replaced with I, which creates an identity matrix of the specified size.
#     The for loop syntax is slightly different, with : used to specify the range and -1 used instead of :-1:1.

#______________________________________________________________________
# PROC missdata                                                        
# PURPOSE: eliminates the rows in y & matrices C, G that correspond to     
#          missing data (NaN) in y                                                                                  
# INPUT    y             vector of observations at time t  (n x 1 )    
#          S             KF system matrices             (structure)
#                        must contain C & G
# OUTPUT   y             vector of observations (reduced)   (# x 1)     
#          C G           KF system matrices     (reduced)   (# x ?)     
#          L             To restore standard dimensions     (n x #)     
#                        where # is the nr of available data in y
#______________________________________________________________________
function MissData(y, C, R)
    ix = findall(.!isnan.(y))  # use .! to broadcast ! to each element of y
    e = eye(length(y))
    L = e[:, ix]

    y = y[ix]
    C = C[ix, :]
    R = R[ix, ix]

    return y, C, R, L
end