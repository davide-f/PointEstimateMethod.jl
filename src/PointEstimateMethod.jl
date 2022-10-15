module PointEstimateMethod

using Distributions
using Combinatorics
using JuMP
using GLPK
using Polynomials

export pem

include("auxiliaries.jl")
include("pem.jl")

end
