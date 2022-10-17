using Revise
using Distributions
using YAML
using PointEstimateMethod


distribution = Normal()  # distribution to consider
N_pem = 3  # number of point estimate models to use

pem_output = pem(distribution, N_pem)


d2 = truncated(Normal(1.0, 0.4), 0.0, +Inf)

pem(d2, 11)