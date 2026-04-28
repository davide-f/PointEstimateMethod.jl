# PointEstimateMethod.jl

**PointEstimateMethod.jl** is a Julia package that implements the well-known Point Estimate Method (PEM) to represent a probability distribution with a discrete set of points that preserve the moments of the original distribution.

## Overview

The Point Estimate Method (PEM) is a technique used in probabilistic analysis to approximate a continuous probability distribution using a finite set of representative points (locations) and their associated probabilities. The selected points preserve the statistical moments (mean, variance, skewness, kurtosis, etc.) of the original distribution up to a desired order.

This package implements the Miller and Rice methodology as described in:

- H.P. Hong, *An efficient point estimate method for probabilistic analysis*, Reliability Engineering and System Safety, 1998, [doi:10.1016/S0951-8320(97)00071-9](https://doi.org/10.1016/S0951-8320(97)00071-9)
- Miller, Allen C., and Thomas R. Rice, *Discrete Approximations of Probability Distributions*, Management Science 29, no. 3 (1983): 352–62, [jstor:2631060](http://www.jstor.org/stable/2631060)

## Quick Start

```julia
using Distributions
using PointEstimateMethod

distribution = Normal()  # distribution to consider
N_pem = 3                # number of point estimate points to use

pem_output = pem(distribution, N_pem)

pem_output.x  # locations of the points
pem_output.p  # probability of each point
```

## Features

- Supports any univariate distribution from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) that provides central moments via `Distributions.moment`.
- Falls back to Monte Carlo sampling for distributions that do not provide a moment function.
- Supports experimental (empirical) distributions provided as vectors of samples.
- Allows user-defined moment functions via the `central_moment_fun` keyword argument.
- Uses [HiGHS](https://highs.dev/) as the default linear programming solver (via [JuMP.jl](https://jump.dev/)).

## Contents

```@contents
Pages = ["index.md", "installation.md", "examples/normal_distribution.md", "examples/general_distribution.md", "examples/multivariate_distribution.md", "api_reference.md"]
Depth = 2
```
