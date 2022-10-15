"""
    pem(d, N; std_fun=std, central_moment_fun=moment)

Point Estimate Method to identify N estimate points for the univariate distribution d.
This function is based on the methodology proposed by:
- 

Parameters
----------
- d :: UnivariateDistribution
    Distribution under interest
- N :: Integer
    Number of desired estimate points
- mean_fun :: Function (optional)
    Function used to calculate the mean of the distribution
- std_fun :: Function (optional)
    Function used to calculate the standard deviation of the distribution
- central_moment_fun :: Function (optional)
    Function used to calculate the central moment of the distribution d

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
        std_fun::Function=Distributions.std,
        central_moment_fun::Function=Distributions.moment,
    )
    
    ## Execution
    ## Solving methodology by https://www.jstor.org/stable/2631060
    
    # std variable
    std_value = std_fun(d)
    
    # lambda i value
    λ = Dict(
        i=>central_moment_fun(d, i)/(std_value^i)
        for i = 1:(2*N)
    )
    λ[0] = 1.0
    
    
    ## 1) Preliminary model to get the coefficients of polynomial described in section 4
    ##    of https://www.jstor.org/stable/2631060
    
    model = Model(GLPK.Optimizer)
    
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
    
    value.(C)
    
    # 2) postprocess the coefficients to obtain the desired locations
    poly_coeffs = [[value(C[i]) for i = 0:N-1]; 1.0]
    
    poly = Polynomial(poly_coeffs)
    
    # get the roots of the polynomial and get the standardized locations of the distribution
    ϵ = roots(poly)
    
    # obtain the probabilities of such locations
    postmodel = Model(GLPK.Optimizer)
    
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
    x = ϵ.*std_value .+ mean_fun(d)
    
    # results = NamedTuple{(:x, :p)}.(zip(x, p))

    return (x=x, p=p)
end