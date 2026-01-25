##################
# Run BM14 Code
##################

cd("/Users/sebastiankrantz/Documents/R/dfms/misc/BM2014/RepFilesBanburaModugnoJAE14")
pwd()

# Data
# using CSV, DataFrames
# # Read CSV file
# data = CSV.read("Data/X_MQ_diff.csv", DataFrame)
# # Print the data
# print(data)
#
# X = Matrix(data[:, 2:end])
# replace!(X, missing => NaN)
# X = Float64.(X)
# X = hcat(X[:, 10:end], X[:, 1:9]) # quarterly variables need to be at the end

# BM14 data as in dfms package
data = CSV.read("../BM14_MQ_12.csv", DataFrame)
X = Matrix(data)
replace!(X, missing => NaN)
X = Float64.(X)

x = Float64.(Matrix(CSV.read("../BM14_MQ_12_imp.csv", DataFrame)))

include("EM_DFM_SS_restrMQ.jl")

# Parameters
P = (r = 2, p = 3, max_iter = 100,
# Building the matrix of restrictions on the loadings for the quarterly variables
    Rconstr = [2 -1 0 0 0;
               3 0 -1 0 0;
               2 0 0 -1 0;
               1 0 0 0 -1],
    q = zeros(4,1),
    nQ = 2 # need to set number of quarterly variables
)


# Run DFM
res = EM_DFM_SS_restrMQ(X, P);

keys(res)

using Plots
plot(res["F"])
