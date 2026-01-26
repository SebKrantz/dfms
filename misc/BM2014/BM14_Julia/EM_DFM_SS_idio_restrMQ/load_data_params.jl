##################
# Run BM14 Code - MQ + AR(1) Idiosyncratic
##################

cd("/Users/sebastiankrantz/Documents/R/dfms/misc/BM2014/BM14_Julia/EM_DFM_SS_idio_restrMQ")
pwd()

# Data
using CSV, DataFrames

# BM14 data as in dfms package
data = CSV.read("../../BM14_MQ_12.csv", DataFrame)
X = Matrix(data)
replace!(X, missing => NaN)
X = Float64.(X)

include("EM_DFM_SS_idio_restrMQ.jl")

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
res = EM_DFM_SS_idio_restrMQ(X, P);

keys(res)

# Check results
println("Number of iterations: ", res["num_iter"])
println("Converged: ", res["converged"])
println("Factor dimensions: ", size(res["F"]))

using Plots
plot(res["F"], title="Estimated Factors (MQ + AR1 Idio)")
