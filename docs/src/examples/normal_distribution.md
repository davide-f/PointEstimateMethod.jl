# Example: Normal Distribution

This example demonstrates the use of `PointEstimateMethod.jl` to approximate a standard Normal distribution using the Point Estimate Method (PEM).

## Setup

```julia
using Distributions
using PointEstimateMethod
```

## Approximating a Standard Normal Distribution

The standard Normal distribution has mean 0 and variance 1. We will use 3 estimate points to represent it.

```julia
distribution = Normal()  # Standard Normal: mean=0, std=1
N_pem = 3                # Number of point estimate points

pem_output = pem(distribution, N_pem)
```

The output `pem_output` is a named tuple with two fields:

- `pem_output.x` — the locations of the representative points
- `pem_output.p` — the probability weights of each point

```julia
println("Locations:     ", pem_output.x)
println("Probabilities: ", pem_output.p)
```

The probabilities should sum to 1 and the weighted mean of the locations should equal the distribution mean (0 for the standard Normal):

```julia
using LinearAlgebra

weighted_mean = dot(pem_output.p, pem_output.x)
println("Weighted mean: ", weighted_mean)  # ≈ 0.0
```

## Using More Points

Increasing the number of estimate points improves the moment-matching accuracy:

```julia
for N in [2, 3, 5, 7]
    out = pem(Normal(), N)
    println("N = $N: locations = $(round.(real.(out.x), digits=4)), probs = $(round.(real.(out.p), digits=4))")
end
```

## Truncated Normal Distribution

The PEM also works with truncated distributions:

```julia
d_trunc = truncated(Normal(1.0, 0.4), 0.0, +Inf)
N_pem = 3

pem_trunc = pem(d_trunc, N_pem)
println("Locations:     ", pem_trunc.x)
println("Probabilities: ", pem_trunc.p)
```

## Custom Solver

By default, PointEstimateMethod.jl uses HiGHS as the LP solver. You can specify a different JuMP-compatible solver using the `optimizer` keyword argument:

```julia
using HiGHS

pem_output = pem(Normal(), 3; optimizer=optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
```
