# Installation

## Requirements

- Julia 1.0 or later

## Installing PointEstimateMethod.jl

PointEstimateMethod.jl is a registered Julia package. You can install it using the Julia package manager from the Julia REPL:

```julia
using Pkg
Pkg.add("PointEstimateMethod")
```

Or by entering `]` in the Julia REPL to open the package manager, then running:

```
pkg> add PointEstimateMethod
```

## Dependencies

PointEstimateMethod.jl automatically installs its dependencies:

- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) — probability distributions
- [JuMP.jl](https://jump.dev/) — mathematical optimization
- [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl) — default LP solver
- [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl) — polynomial root finding
- [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) — combinatorics utilities

## Verifying the Installation

After installation, you can verify that the package works correctly by running:

```julia
using Distributions
using PointEstimateMethod

pem_output = pem(Normal(), 3)
println("Locations: ", pem_output.x)
println("Probabilities: ", pem_output.p)
```
