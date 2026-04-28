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
        d::Vector{<:UnivariateDistribution},
        N::Integer;
        mean_fun::Function=Distributions.mean,
        central_moment_fun::Function=Distributions.moment,
        montecarlo_sampling::Integer= 1000000,
        optimizer=DEFAULT_SOLVER,
    )

    K = length(d)
    
    if all(hasmethod(central_moment_fun, Tuple{typeof(d[k]), Int}) for k in 1:K)
    
        # central moments
        m_list = Dict(
            (k,i)=>central_moment_fun(d[k], i)
            for k = 1:K for i = 1:(2*N)
        )

        mean_values = [mean_fun(d[k]) for k in 1:K]

        return pem(mean_values, m_list, N; optimizer=optimizer)
    else
        @info """Function $(string(central_moment_fun)) does not have a direct implementation for Distribution $(string(d)). Perform Monte Carlo sempling over the distribution with $montecarlo_sampling points"""
        sampled_set = zeros(K, montecarlo_sampling)
        for k = 1:K
            sampled_set[k, :] = rand(d[k], montecarlo_sampling)
        end

        return pem(sampled_set, N; optimizer=optimizer)
    end
end


"""
    pem(d, N; mean_fun=mean, central_moment_fun=moment, optimizer=HiGHS.Optimizer)

Point Estimate Method to identify N estimate points for an experimental multi-variate distribution
represented by the array of elements d.
The first dimension of the array is the number of variables, while the second dimension is the number of samples for each variable.

Parameters
----------
- d :: Array{<:Real, 2}
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
    - x :: Array
        Return the location points
    - p :: Array
        Return the probability of each point

"""
function pem(
        d::Array{<:Real, 2},
        N::Integer;
        mean_fun::Function=Distributions.mean,
        central_moment_fun::Function=Distributions.moment,
        optimizer=DEFAULT_SOLVER,
    )
    
    ## Execution
    ## Solving methodology by https://www.jstor.org/stable/2631060

    K = size(d, 1)

    mean_values = [mean_fun(d[k,:]) for k in 1:K]
    
    # moments
    m_list = Dict(
        (k, i)=>central_moment_fun(d[k,:], i)
        for k = 1:K for i = 1:(2*N)
    )

    return pem(mean_values, m_list, N; optimizer=optimizer)
end


"""
    pem(mean_values, d, m_list, N; optimizer=HiGHS.Optimizer)

Point Estimate Method to identify estimate points for K independent univariate distributions with moments given by m_list and mean values given by mean_values.

This function is based on the methodology proposed by:
- H.P.Hong, An efficient point estimate method for probabilistic analysis, Reliability Engineering and System Safety, 1998, https://doi.org/10.1016/S0951-8320(97)00071-9
- Miller, Allen C., and Thomas R. Rice. “Discrete Approximations of Probability Distributions.” Management Science 29, no. 3 (1983): 352–62. http://www.jstor.org/stable/2631060.

Parameters
----------
- mean_value :: Vector{<:Real}
    Mean value of the distribution
- m_list :: Dict{Tuple{<:Integer,<:Integer}, <:Real}
    Dictionary representing the central moments of the distribution.
    The keys of the dictionary are tuples (k, m) where k is the distribution index and m is the moment order.
    The value corresponds to the value of the moment.
    Note: they are central moments referred to the mean. As such, the moment of order 1 is 0.0.
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
        mean_values::Vector{<:Real},
        m_list::Dict{<:Tuple{<:Integer,<:Integer}, <:Real},
        N::Integer;
        optimizer=DEFAULT_SOLVER,
    )

    # number of distributions
    K = length(mean_values)
    
    # ensure consistency of the input
    expected_keys = Set((k, i) for k in 1:K for i in 1:(2*N))
    @assert Set(keys(m_list)) == expected_keys "The input moment dictionary does not match the expected index in th form (k,m) where k is the distribution index and m is the moment order."
    
    ## Execution
    ## Solving methodology by https://www.jstor.org/stable/2631060
    
    # lambda i value
    # The moment of order 0 has value 1/K, as we are considering K independent distributions
    λ = Dict(
        (k,i)=>((i==0) ? 1. / K : m_list[(k,i)])
        for k=1:K for i = 0:2*N
    )
    
    ## 1) Preliminary model to get the coefficients of polynomial described in section 4
    ##    of https://www.jstor.org/stable/2631060
    
    model = Model(optimizer)
    
    # coefficients of auxiliary polynomial \sum_{k=0}^N C_k x^k = π(x) = (x - x_1) ... (x - x_N) for each distribution
    @variable(model, C[k=1:K,i=0:N-1])

    @constraint(
        model,
        aux_poly_balance[k=1:K,i=0:N-1],
        sum(
            C[k,p]*λ[k,p+i]
            for p=0:N-1
        ) == -λ[k,N+i]
    )
    
    # Determine coefficients
    optimize!(model)
    
    # 2) postprocess the coefficients to obtain the desired locations
    ϵ = zeros(K, N)

    for k in 1:K
        # get the coefficients of the polynomial
        poly_coeffs = [[value(C[k,i]) for i = 0:N-1]; 1.0]
        # get the roots of the polynomial and get the standardized locations of the distribution
        poly = Polynomials.Polynomial(poly_coeffs)
        ϵ[k,:] = Polynomials.roots(poly)
    end
    
    # obtain the probabilities of such locations
    postmodel = Model(optimizer)
    
    @variable(postmodel, probabilities[k=1:K,j=1:N])

    # balance the moments for each distribution
    @constraint(
        postmodel,
        balance_moments[k=1:K,i=1:(2*N-1)],
        sum(probabilities[k,j] * ϵ[k,j]^i for j in 1:N) == λ[k,i]
    )
    
    # set the probabilities to sum to 1/K for each distribution, as we are considering K independent distributions
    @constraint(
        postmodel,
        balance_probabilities[k=1:K],
        sum(probabilities[k,:]) == 1. / K
    )

    # optimize the model
    optimize!(postmodel)
    
    # get probabilities
    p = reshape(value.(probabilities), K, N)

    # get locations of the estimated points
    x = ϵ .+ reshape(mean_values, K, 1)
    
    # results = NamedTuple{(:x, :p)}.(zip(x, p))
    return (x=x, p=p)
end
