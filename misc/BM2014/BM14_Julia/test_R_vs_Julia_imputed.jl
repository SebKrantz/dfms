###########################################
# Test Julia vs R DFM Implementations
# Using pre-imputed data to eliminate preprocessing differences
###########################################

cd(@__DIR__)
pwd()

using CSV, DataFrames
using LinearAlgebra
using Statistics

###########################################
# Load pre-imputed data (from R's tsnarmimp)
###########################################
println("=== Loading Pre-imputed Data ===")

data_imp = CSV.read("../BM14_MQ_12_imp.csv", DataFrame)
X_imp = Matrix(data_imp)
X_imp = Float64.(X_imp)

data_raw = CSV.read("../BM14_MQ_12.csv", DataFrame)
X_raw = Matrix(data_raw)
replace!(X_raw, missing => NaN)
X_raw = Float64.(X_raw)

println("Imputed data dimensions: ", size(X_imp))
println("Raw data dimensions: ", size(X_raw))
println("NaN count in imputed: ", sum(isnan.(X_imp)))
println("NaN count in raw: ", sum(isnan.(X_raw)))

###########################################
# Compare initial PCA on imputed data
###########################################
println("\n=== Initial PCA Comparison ===")

# Standardize imputed data
Mx = mean(X_imp, dims=1)
Wx = std(X_imp, dims=1)
X_std = (X_imp .- Mx) ./ Wx

# Compute eigendecomposition
Sigma = cov(X_std)
eig = eigen(Sigma)

# Get largest 2 eigenvalues (last 2 in ascending order)
r = 2
d = eig.values[end-r+1:end]
v = eig.vectors[:, end-r+1:end]

println("Largest 2 eigenvalues: ", reverse(d))
println("\nFirst eigenvector (corresponding to largest eigenvalue):")
println("  v[:, end] = ", round.(v[:, end], digits=6))
println("\nSecond eigenvector:")
println("  v[:, end-1] = ", round.(v[:, end-1], digits=6))

###########################################
# Compare with R's output
###########################################
println("\n=== Expected from R ===")
println("R C matrix (first 5 rows):")
println("  ip_total:                  0.673461,  0.089035")
println("  ip_en_2:                   0.000628,  0.051399")
println("  ret_turnover_defl:         0.167502,  0.098082")
println("  ecs_cons_sit_over_next_12: 0.129640,  0.049497")
println("  ecs_ret_tr_exp_bus:        0.028561, -0.069781")

# For MQ model, C = [v, zeros(N, r*(pC-1))] initially
# where pC = max(p, rC) = max(3, 5) = 5
println("\n=== Julia Initial C (from eigenvectors) ===")
println("v (first 5 rows):")
for i in 1:5
    println("  ", names(data_imp)[i], ": ", round(v[i, end], digits=6), ", ", round(v[i, end-1], digits=6))
end

###########################################
# Check R's PCA implementation
###########################################
println("\n=== Check: R uses prcomp which centers and scales ===")
# R's prcomp returns loadings such that X %*% loadings = scores
# Let's verify our eigenvectors match the expected pattern

# PCA via SVD (like R's prcomp)
F = svd(X_std)
loadings_svd = F.V[:, 1:r]  # First r right singular vectors

println("SVD loadings (first 5 rows, should match R):")
for i in 1:5
    println("  ", names(data_imp)[i], ": ", round(loadings_svd[i, 1], digits=6), ", ", round(loadings_svd[i, 2], digits=6))
end

# Note: sign may be flipped - let's check correlation with R's values
# R C[1,1] = 0.673461, our loading should be similar in magnitude
println("\nMagnitude check:")
println("  |loadings_svd[1,1]| = ", abs(loadings_svd[1, 1]))
println("  R |C[1,1]| = 0.673461")
