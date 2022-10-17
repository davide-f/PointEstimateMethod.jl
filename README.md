![CI](https://github.com/davide-f/PointEstimateMethod.jl/actions/workflows/CI.yml/badge.svg)


# PointEstimateMethod.jl

This repository implements the well known Point Estimate Method (PEM) to represent a distribution
with a number of points that preserve the moments of the original distribution.

The methodology implements the Miller and Rice methodology as described at:
- H.P.Hong, An efficient point estimate method for probabilistic analysis, Reliability Engineering and System Safety, 1998, https://doi.org/10.1016/S0951-8320(97)00071-9
- Miller, Allen C., and Thomas R. Rice. “Discrete Approximations of Probability Distributions.” Management Science 29, no. 3 (1983): 352–62. http://www.jstor.org/stable/2631060.

Example of usage is as follow:

```julia

using Distributions
using PointEstimateMethod

distribution = Normal()  # distribution to consider
N_pem = 3  # number of point estimate models to use

pem_output = pem(distribution, N_pem)

pem_output.x  # locations of the points
pem_output.p  # probability of each point
```