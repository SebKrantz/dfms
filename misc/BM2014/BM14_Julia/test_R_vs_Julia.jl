###########################################
# Test Julia vs R DFM Implementations
# Compare EM_DFM_SS_restrMQ and News_DFM_ML
###########################################

cd(@__DIR__)
pwd()

using CSV, DataFrames
using LinearAlgebra
using Statistics
using NaNStatistics

###########################################
# Test 1: MQ Model (EM_DFM_SS_restrMQ)
###########################################
println("=== Test 1: MQ Model (r=2, p=3) ===")

include("EM_DFM_SS_restrMQ/EM_DFM_SS_restrMQ.jl")

# Load BM14 data (same as used in R)
data = CSV.read("../BM14_MQ_12.csv", DataFrame)
X = Matrix(data)
replace!(X, missing => NaN)
X = Float64.(X)

println("Data dimensions: ", size(X))
println("Column names: ", names(data))

# Parameters matching R
P_mq = (r = 2, p = 3, max_iter = 100,
    Rconstr = [2 -1 0 0 0;
               3 0 -1 0 0;
               2 0 0 -1 0;
               1 0 0 0 -1],
    q = zeros(4,1),
    nQ = 2  # invest, productivity are quarterly
)

# Run MQ model
res_mq = EM_DFM_SS_restrMQ(X, P_mq)

println("Converged: ", res_mq["converged"])
println("Iterations: ", res_mq["num_iter"])
println("Factor dimensions: ", size(res_mq["F"]))

println("\n--- Parameters comparison ---")
println("C matrix dimensions: ", size(res_mq["C"]))
println("A matrix dimensions: ", size(res_mq["A"]))
println("Q matrix dimensions: ", size(res_mq["Q"]))

println("\nC matrix (first 5 rows, first 2 cols):")
display(round.(res_mq["C"][1:5, 1:2], digits=6))

###########################################
# Test 2: News Decomposition (MQ model)
###########################################
println("\n\n=== Test 2: News Decomposition (MQ) ===")

# Create old and new vintages
X_new = copy(X)
X_old = copy(X)

# Simulate data releases: remove recent monthly observations (same as R)
X_old[350, 1] = NaN  # ip_total
X_old[351, 1] = NaN

println("Releases at t=350,351 for variable 1 (ip_total)")

# Fit models on both vintages
res_old = EM_DFM_SS_restrMQ(X_old, P_mq)
res_new = EM_DFM_SS_restrMQ(X_new, P_mq)

println("Old model converged: ", res_old["converged"], " (", res_old["num_iter"], " iterations)")
println("New model converged: ", res_new["converged"], " (", res_new["num_iter"], " iterations)")

# Set up news decomposition
include("News/News_DFM_ML.jl")

# Create Q dict for news function (using OLD model parameters, like MATLAB/R)
Q_news = Dict(
    "C" => res_old["C"],
    "R" => res_old["R"],
    "A" => res_old["A"],
    "Q" => res_old["Q"],
    "Z_0" => res_old["Z_0"],
    "V_0" => Matrix(res_old["V_0"]),
    "Mx" => res_old["Mx"],
    "Wx" => res_old["Wx"],
    "Groups" => vcat(repeat(["M"], 10), ["Q", "Q"]),  # 10 Monthly + 2 Quarterly
    "Series" => string.(names(data))
)

# Run news decomposition
t_fcst = 355
v_news = 3  # ret_turnover_defl (variable 3, matching R)

println("\nTarget time: ", t_fcst)
println("Target variable: ", v_news, " (", names(data)[v_news], ")")

y_old, y_new, groupnews, singlenews, gain, gainSer, actual, fore, filt =
    News_DFM_ML(X_old, X_new, Q_news, t_fcst, v_news)

println("\n--- News Results (Julia) ---")
println("Old forecast: ", y_old)
println("New forecast: ", y_new)
println("Revision: ", y_new - y_old)
println("Sum of news: ", sum(singlenews))
println("Difference (should be ~0): ", (y_new - y_old) - sum(singlenews))

println("\nNon-zero news contributions:")
nz = findall(singlenews .!= 0)
for i in nz
    println("  Variable ", i, " (", names(data)[i], "): ", singlenews[i])
end

###########################################
# Test 3: MQ + AR1 Model
###########################################
println("\n\n=== Test 3: MQ + AR1 Idiosyncratic Model ===")

include("EM_DFM_SS_idio_restrMQ/EM_DFM_SS_idio_restrMQ.jl")

# Parameters for MQ + AR1
P_mq_ar1 = (r = 2, p = 3, max_iter = 100,
    Rconstr = [2 -1 0 0 0;
               3 0 -1 0 0;
               2 0 0 -1 0;
               1 0 0 0 -1],
    q = zeros(4,1),
    nQ = 2
)

# Run MQ + AR1 model
try
    res_mq_ar1 = EM_DFM_SS_idio_restrMQ(X, P_mq_ar1)

    println("Converged: ", res_mq_ar1["converged"])
    println("Iterations: ", res_mq_ar1["num_iter"])

    println("\nC matrix (first 5 rows, first 2 cols):")
    display(round.(res_mq_ar1["C"][1:5, 1:2], digits=6))

    # Extract AR1 coefficients from A matrix
    # For MQ+AR1, A is block diagonal: [A_factors | B_monthly | B_quarterly]
    r = P_mq_ar1.r
    pC = max(P_mq_ar1.p, 5)
    rpC = r * pC
    nM = size(X, 2) - P_mq_ar1.nQ

    A_full = res_mq_ar1["A"]
    B_monthly = diag(A_full[rpC+1:rpC+nM, rpC+1:rpC+nM])

    println("\nAR1 coefficients for monthly variables:")
    for i in 1:nM
        println("  ", names(data)[i], ": ", round(B_monthly[i], digits=4))
    end

catch e
    println("Error fitting MQ+AR1 model: ", e)
    println(stacktrace(catch_backtrace()))
end

###########################################
# Summary Comparison with R
###########################################
println("\n\n=== Summary: Comparison with R ===")
println("R results (from test_R_vs_Julia.R):")
println("  MQ model: C[1,1] = 0.673461, C[1,2] = 0.089035")
println("  News: y_old = 0.002362266, y_new = -0.002901152")
println("  News: revision = -0.005263418, sum(singlenews) = -0.005263418")
println("\nJulia results:")
println("  MQ model: C[1,1] = ", round(res_mq["C"][1,1], digits=6), ", C[1,2] = ", round(res_mq["C"][1,2], digits=6))
println("  News: y_old = ", y_old, ", y_new = ", y_new)
println("  News: revision = ", y_new - y_old, ", sum(singlenews) = ", sum(singlenews))
