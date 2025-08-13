# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 Davide Fioriti
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8

"""
    pem(d, N; mean_fun=mean, central_moment_fun=moment, optimizer=HiGHS.Optimizer)

Point Estimate Method to identify N estimate points for the univariate distribution d.

Parameters
----------
- d :: UnivariateDistribution
    Distribution under interest
- N :: Integer
    Number of desired estimate points
- mean_fun :: Function (optional)
    Function used to calculate the mean value of the distribution
- central_moment_fun :: Function (optional)
    Function used to calculate the central moment of the distribution d
- montecarlo_sampling :: Integer (optional, default 1e6)
    Number of Monte Carlo samples used in the sampling procedure if a non-specific moment function is available
- optimizer (optional)
    JuMP optimizer for executing the optimization

Returns
-------
- (x, p) :: NamedTuple
    - x :: Vector
        Return the location points
    - p :: Vector
        Return the probability of each point

"""
function pem(
        d::UnivariateDistribution,
        N::Integer;
        mean_fun::Function=Distributions.mean,
        central_moment_fun::Function=Distributions.moment,
        montecarlo_sampling::Integer= 1000000,
        optimizer=DEFAULT_SOLVER,
    )
    
    if hasmethod(central_moment_fun, Tuple{typeof(d), Int})
    
        # lambda i value
        m_list = Dict(
            i=>central_moment_fun(d, i)
            for i = 1:(2*N)
        )

        return pem(mean_fun(d), m_list, N; optimizer=optimizer)
    else
        @info """Function $(string(central_moment_fun)) does not have a direct implementation for Distribution $(string(d)). Perform Monte Carlo sempling over the distribution with $montecarlo_sampling points"""
        sampled_set = rand(d, montecarlo_sampling)

        return pem(sampled_set, N; optimizer=optimizer)
    end
end

"""
    pem(d, N; mean_fun=mean, central_moment_fun=moment, optimizer=HiGHS.Optimizer)

Point Estimate Method to identify N estimate points for an experimental distribution
represented by the vector of elements d.

Parameters
----------
- d :: Vector
    Distribution under interest
- N :: Integer
    Number of desired estimate points
- mean_fun :: Function (optional)
    Function used to calculate the mean value of the distribution
- central_moment_fun :: Function (optional)
    Function used to calculate the central moment of the distribution d
- optimizer (optional)
    JuMP optimizer for executing the optimization

Returns
-------
- (x, p) :: NamedTuple
    - x :: Vector
        Return the location points
    - p :: Vector
        Return the probability of each point

"""
function pem(
        d::Vector,
        N::Integer;
        mean_fun::Function=Distributions.mean,
        central_moment_fun::Function=Distributions.moment,
        optimizer=DEFAULT_SOLVER,
    )
    
    ## Execution
    ## Solving methodology by https://www.jstor.org/stable/2631060
    
    # lambda i value
    m_list = Dict(
        i=>central_moment_fun(d, i)
        for i = 1:(2*N)
    )

    return pem(mean_fun(d), m_list, N; optimizer=optimizer)
end


"""
    pem(mean_value, d, m_list, N; optimizer=HiGHS.Optimizer)

Point Estimate Method to identify N estimate points for the univariate distribution d.
This function is based on the methodology proposed by:
- H.P.Hong, An efficient point estimate method for probabilistic analysis, Reliability Engineering and System Safety, 1998, https://doi.org/10.1016/S0951-8320(97)00071-9
- Miller, Allen C., and Thomas R. Rice. “Discrete Approximations of Probability Distributions.” Management Science 29, no. 3 (1983): 352–62. http://www.jstor.org/stable/2631060.

Parameters
----------
- mean_value
    Mean value of the distribution
- m_list :: Dict
    Dictionary representing the moments of the distribution.
    The keys of the dictionary shall go from 0 to N and the value corresponds to the value
    of the moment
- N :: Integer
    Number of desired estimate points
- optimizer (optional)
    JuMP optimizer for executing the optimization

Returns
-------
- (x, p) :: NamedTuple
    - x :: Vector
        Return the location points
    - p :: Vector
        Return the probability of each point

"""
function pem(
        mean_value,
        m_list::Dict,
        N::Integer;
        optimizer=DEFAULT_SOLVER,
    )
    
    @assert Set(keys(m_list)) == Set(1:2*N) "The input moment dictionary does not match the expected index 1,...,2*N"
    
    ## Execution
    ## Solving methodology by https://www.jstor.org/stable/2631060
    
    # lambda i value
    λ = Dict(i=>m_list[i] for i = 1:2*N)
    λ[0] = 1.0    
    
    ## 1) Preliminary model to get the coefficients of polynomial described in section 4
    ##    of https://www.jstor.org/stable/2631060
    
    model = Model(optimizer)
    
    # coefficients of auxiliary polynomial \sum_{k=0}^N C_k x^k = π(x) = (x - x_1) ... (x - x_N)
    @variable(model, C[i=0:N-1])
    
    @constraint(
        model,
        aux_poly_balance[i=0:N-1],
        sum(
            C[p]*λ[p+i]
            for p=0:N-1
        ) == -λ[N+i]
    )
    
    # Determine coefficients
    optimize!(model)
    
    # 2) postprocess the coefficients to obtain the desired locations
    poly_coeffs = [[value(C[i]) for i = 0:N-1]; 1.0]
    
    poly = Polynomials.Polynomial(poly_coeffs)
    
    # get the roots of the polynomial and get the standardized locations of the distribution
    ϵ = Polynomials.roots(poly)
    
    # obtain the probabilities of such locations
    postmodel = Model(optimizer)
    
    @variable(postmodel, 0 <= probabilities[j=1:N] <= 1.0)
    
    @constraint(
        postmodel,
        balance_moments[i=0:(2*N-1)],
        sum(probabilities[j] * ϵ[j]^i for j in 1:N) == λ[i]
    )
    
    optimize!(postmodel)
    
    # get probabilities
    p = value.(probabilities)

    # get locations of the estimated points
    x = ϵ .+ mean_value
    
    # results = NamedTuple{(:x, :p)}.(zip(x, p))
    return (x=x, p=p)
end