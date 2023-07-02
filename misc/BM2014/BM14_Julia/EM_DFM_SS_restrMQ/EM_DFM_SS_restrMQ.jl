###########################################
# Julia Implementation of Mixed Frequency
# DFM with WN Observation Residuals
###########################################

# Installation and loading
# using Pkg
# Pkg.add("NaNStatistics") # nanmean etc.
# Pkg.add("DSP") # filt()
# Pkg.add("DataInterpolations") # CubicSpline()
# Pkg.add("BlockDiagonals") # BlockDiagonal()
using Statistics
using NaNStatistics
using DSP
using DataInterpolations
using LinearAlgebra
using BlockDiagonals

eye(n) = Matrix{Float64}(I, n, n)


mutable struct Nanopt
    method
    k
end

# Functions 
include("runKF.jl")
include("remNaNs_spline.jl")
include("Procedures.jl")

# Parameters
P = (r = 4, p = 2, max_iter = 100,
# Building the matrix of restrictions on the loadings for the quarterly variables
    Rconstr = [2 -1 0 0 0;
               3 0 -1 0 0;
               2 0 0 -1 0;
               1 0 0 0 -1],
    q = zeros(4,1),
    nQ = 2 # need to set number of quarterly variables
)


function EM_DFM_SS_restrMQ(X, P)

    # X: T x N panel of data
    # P: structure containing settings of the model and estimation 

    thresh = 1e-4
    r = P.r # number of factors
    p = P.p # number of lags in the factor VAR
    max_iter = P.max_iter # maximal number of EM iterations
    R_mat = P.Rconstr # matrix of restrictions on the loadings for the quarterly variables
    q = P.q
    nQ = P.nQ # number of quarterly variables

    #--------------------------------------------------------------------------
    # Preparation of the data
    #--------------------------------------------------------------------------
    T, N = size(X)

    # Standardise x
    Mx = nanmean(X, dims = 1)
    Wx = nanstd(X, dims = 1)

    # Standardize X using Mx and Wx  
    xNaN = (X .- Mx) ./ Wx
    # xNaN = X;

    #--------------------------------------------------------------------------
    # Initial Conditions
    #--------------------------------------------------------------------------

    #Removing missing values (for initial estimators)

    optNaN = Nanopt(2, 3) # Remove leading and closing zeros

    A, C, Q, R, Z_0, V_0 = InitCond(xNaN, r, p, optNaN, R_mat, q, nQ)

    # some auxiliary variables for the iterations
    previous_loglik = -Inf
    num_iter = 0
    LL = -Inf
    converged = false

    # y for the estimation is WITH missing data
    y = transpose(xNaN)

    #--------------------------------------------------------------------------
    #THE EM LOOP
    #--------------------------------------------------------------------------

    #The model can be written as
    #y = C*Z + e;
    #Z = A*Z(-1) + v
    #where y is NxT, Z is (pr)xT, etc

    #remove the leading and ending nans for the estimation
    optNaN.method = 3
    y_est = remNaNs_spline(xNaN, optNaN)[1]'

    decrease = zeros(max_iter+1)
    LL = zeros(max_iter+1)

    while num_iter < max_iter && !converged 
        C_new, R_new, A_new, Q_new, loglik, Z_0, V_0 = EMstep(y_est, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ) 
        
        C = C_new
        R = R_new
        A = A_new
        Q = Q_new

        if num_iter > 2
            converged, decrease[num_iter+1] = em_converged(loglik, previous_loglik, thresh, true)
        end
        
        LL[num_iter+1] = loglik
        previous_loglik = loglik
        num_iter += 1
    end

    #final run of the Kalman filter
    Zsmooth = runKF(y, A, C, Q, R, Z_0, V_0)[1]'
    Res = Dict()
    Res["x_sm"] = Zsmooth[2:end,:] * C'
    Res["X_sm"] = Wx .* Res["x_sm"] .+ Mx
    Res["F"] = Zsmooth[2:end, 1:r]

    #--------------------------------------------------------------------------
    #   Loading the structure with the results
    #--------------------------------------------------------------------------
    Res["C"] = C
    Res["R"] = R
    Res["A"] = A
    Res["Q"] = Q
    Res["Z_0"] = Z_0
    Res["V_0"] = V_0
    Res["r"] = r
    Res["p"] = p
    Res["Mx"] = Mx
    Res["Wx"] = Wx
    Res["loglik"]   = LL
    Res["num_iter"] = num_iter
    Res["converged"] = converged;
    Res["decreases"] = any(decrease .== 1)
    Res["decr_min_iter"] = minimum(findall(decrease .== 1))

    return Res
end
