# # Example: Multivariate Distribution
#
# This example demonstrates the use of `PointEstimateMethod.jl` to approximate a set of independent univariate distributions jointly using the multivariate PEM.
#
# ## Setup

using Distributions
using PointEstimateMethod

# ## Approximating Two Independent Distributions
#
# We consider two independent distributions:
# - A Normal distribution: ``X_1 \sim \mathcal{N}(2.0, 0.5)``
# - A LogNormal distribution: ``X_2 \sim \text{LogNormal}(0.5, 0.3)``

d1 = Normal(2.0, 0.5)
d2 = LogNormal(0.5, 0.3)

distributions = [d1, d2]
N_pem = 3  # Number of point estimate points per distribution

pem_output = pem(distributions, N_pem)

# The output `pem_output` is a named tuple with two fields:
#
# - `pem_output.x` — a `K × N` matrix of location points (one row per distribution)
# - `pem_output.p` — a `K × N` matrix of probability weights (one row per distribution)

println("Locations matrix (K × N):")
println(pem_output.x)

println("Probabilities matrix (K × N):")
println(pem_output.p)

# ## Verifying Moment Preservation
#
# The sum of probabilities for each distribution should equal `1/K`:

K = length(distributions)
for k in 1:K
    println("Sum of probabilities for distribution $k: ", sum(pem_output.p[k, :]))  # ≈ 1/K
end

# The weighted mean for each distribution should match the true mean:

for k in 1:K
    weighted_mean = sum(pem_output.p[k, j] * pem_output.x[k, j] for j in 1:N_pem) / sum(pem_output.p[k, :])
    println("Distribution $k — PEM weighted mean: $(real(weighted_mean)), true mean: $(mean(distributions[k]))")
end

# ## Using Sample Data
#
# The multivariate PEM also accepts a `K × M` matrix of samples, where `K` is the number of variables and `M` is the number of samples:

using Random
Random.seed!(42)

M = 10000
samples = vcat(
    rand(Normal(2.0, 0.5), 1, M),
    rand(LogNormal(0.5, 0.3), 1, M)
)  # 2 × 10000 matrix

pem_samples = pem(samples, 3)

println("Locations from samples:")
println(pem_samples.x)
println("Probabilities from samples:")
println(pem_samples.p)
