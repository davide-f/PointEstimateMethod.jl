# Example: General (Empirical) Distribution

This example demonstrates the use of `PointEstimateMethod.jl` to approximate an arbitrary experimental distribution represented as a vector of samples.

## Setup

```julia
using Distributions
using PointEstimateMethod
using Random
using Statistics
```

## Generating Sample Data

We start by generating a sample dataset from a known distribution. In practice, this data could come from measurements or simulations.

```julia
Random.seed!(42)

# Generate samples from a LogNormal distribution
samples = rand(LogNormal(0.5, 0.3), 10000)
```

## Applying the Point Estimate Method

When the input is a `Vector`, PEM automatically computes the sample moments and finds the representative points:

```julia
N_pem = 3  # Number of point estimate points

pem_output = pem(samples, N_pem)

println("Locations:     ", pem_output.x)
println("Probabilities: ", pem_output.p)
```

The output is the same named tuple structure as for parametric distributions:

- `pem_output.x` — the representative point locations
- `pem_output.p` — the probability weights (should sum to approximately 1)

## Verifying Moment Preservation

We can verify that the PEM output correctly preserves the mean of the sample distribution:

```julia
using LinearAlgebra

sample_mean = mean(samples)
pem_mean = dot(pem_output.p, pem_output.x)

println("Sample mean: ", round(sample_mean, digits=4))
println("PEM mean:    ", round(real(pem_mean), digits=4))  # should be close to sample_mean
```

## Using Custom Moment Functions

You can provide custom functions to compute the mean and central moments of the distribution
using the keyword arguments `mean_fun` and `central_moment_fun`:

```julia
using Statistics

custom_mean(d, args...) = mean(d)
custom_moment(d, k) = mean((d .- mean(d)) .^ k)

pem_custom = pem(samples, 3; mean_fun=custom_mean, central_moment_fun=custom_moment)
println("Custom PEM locations:     ", pem_custom.x)
println("Custom PEM probabilities: ", pem_custom.p)
```

## Providing Moments Directly

If you already know the mean and central moments of the distribution, you can provide them directly using the three-argument form of `pem`:

```julia
mean_value = mean(samples)
m_list = Dict(i => mean((samples .- mean_value) .^ i) for i in 1:6)

pem_direct = pem(mean_value, m_list, 3)
println("Direct PEM locations:     ", pem_direct.x)
println("Direct PEM probabilities: ", pem_direct.p)
```

This is useful when moments are obtained from analytical expressions or external sources.
