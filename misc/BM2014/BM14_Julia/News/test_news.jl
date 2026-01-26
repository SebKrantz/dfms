##################
# Test News Decomposition
##################

cd("/Users/sebastiankrantz/Documents/R/dfms/misc/BM2014/BM14_Julia/News")
pwd()

using CSV, DataFrames
using LinearAlgebra
using Statistics
using NaNStatistics

# First run the MQ model to get parameters
include("../EM_DFM_SS_restrMQ/EM_DFM_SS_restrMQ.jl")

# BM14 data as in dfms package
data = CSV.read("../../BM14_MQ_12.csv", DataFrame)
X = Matrix(data)
replace!(X, missing => NaN)
X = Float64.(X)

# Parameters
P = (r = 2, p = 3, max_iter = 100,
    Rconstr = [2 -1 0 0 0;
               3 0 -1 0 0;
               2 0 0 -1 0;
               1 0 0 0 -1],
    q = zeros(4,1),
    nQ = 2
)

# Create old and new vintages
X_new = copy(X)
X_old = copy(X)

# Simulate a data release: set some recent values to NaN in old vintage
X_old[350, 1] = NaN  # Remove a monthly observation
X_old[351, 1] = NaN

# Run model on new data
res_new = EM_DFM_SS_restrMQ(X_new, P)

# Create Q dict for news function
Q = Dict(
    "C" => res_new["C"],
    "R" => res_new["R"],
    "A" => res_new["A"],
    "Q" => res_new["Q"],
    "Z_0" => res_new["Z_0"],
    "V_0" => Matrix(res_new["V_0"]),
    "Mx" => res_new["Mx"],
    "Wx" => res_new["Wx"],
    "Groups" => repeat(["M"], size(X, 2) - 2) |> x -> vcat(x, ["Q", "Q"]),  # Monthly + 2 Quarterly
    "Series" => ["v$i" for i in 1:size(X, 2)]
)

# Run news decomposition
include("News_DFM_ML.jl")

t_fcst = 355  # Target forecast time
v_news = 3    # Target variable (a monthly variable)

y_old, y_new, groupnews, singlenews, gain, gainSer, actual, fore, filt =
    News_DFM_ML(X_old, X_new, Q, t_fcst, v_news)

println("Old forecast: ", y_old)
println("New forecast: ", y_new)
println("Revision: ", y_new - y_old)
println("Sum of news: ", sum(singlenews))
println("Difference (should be ~0): ", (y_new - y_old) - sum(singlenews))

println("\nNon-zero news contributions:")
for i in findall(singlenews .!= 0)
    println("  Variable $i: ", singlenews[i])
end
