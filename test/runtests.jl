using Test
using Distributions
using YAML
using PointEstimateMethod

include("Examples.jl")

BASE_FOLDER = dirname(dirname(pathof(PointEstimateMethod)))

ATOL_TEST = 1e-4
RTOL_TEST = 1e-4
SEED = 0

# function to test the approximation of vectors
function approx_test(x_test::Vector, x_validation::Vector, rtol=RTOL_TEST, atol=ATOL_TEST)
    length(x_test) != length(x_validation) && return false
    return all(isapprox(xt, xv; rtol=rtol, atol=atol) for (xt, xv) in zip(x_test, x_validation))
end

# function to test the approximation of matrices
function approx_test(x_test::Matrix, x_validation::Matrix, rtol=RTOL_TEST, atol=ATOL_TEST)
    size(x_test) != size(x_validation) && return false
    return all(isapprox(xt, xv; rtol=rtol, atol=atol) for (xt, xv) in zip(x_test, x_validation))
end

"Function to test examples"
function test_example(example_name, testing_function, args...)

    # for reproducibility
    Distributions.Random.seed!(SEED)
    
    # calculate simulations
    calc_solution = testing_function(args...)

    path_solution = joinpath(BASE_FOLDER, "test", "testcases", string(testing_function), example_name * ".yml")

    
    if calc_solution.x isa Vector
        # Univariate distribution

        I_calc = sortperm(calc_solution.x)
        x_sorted = calc_solution.x[I_calc]
        p_sorted = calc_solution.p[I_calc]

    else
        # Multivariate distribution
        K = size(calc_solution.x, 1)
        I_calc = [sortperm(calc_solution.x[k, :]) for k in 1:K]
        x_sorted = Matrix(reduce(hcat, [calc_solution.x[k, I_calc[k]] for k in 1:K])')
        p_sorted = Matrix(reduce(hcat, [calc_solution.p[k, I_calc[k]] for k in 1:K])')
    end
    
    if isfile(path_solution)
        # if the file exists run tests
        proven_solution = YAML.load_file(path_solution)

        x_validation = proven_solution["x"] isa Vector{<:Vector} ? Matrix(reduce(hcat, proven_solution["x"])') : proven_solution["x"]
        p_validation = proven_solution["p"] isa Vector{<:Vector} ? Matrix(reduce(hcat, proven_solution["p"])') : proven_solution["p"]
        
        @test approx_test(x_sorted, x_validation)
        @test approx_test(p_sorted, p_validation)
    else
        # otherwise create the tests
        mkpath(dirname(path_solution))

        x_for_yaml = x_sorted isa Matrix ? [x_sorted[i, :] for i in axes(x_sorted, 1)] : x_sorted
        p_for_yaml = p_sorted isa Matrix ? [p_sorted[i, :] for i in axes(p_sorted, 1)] : p_sorted


        dict_calc_solution = Dict(
            "x"=>x_for_yaml,
            "p"=>p_for_yaml,
        )

        YAML.write_file(path_solution, dict_calc_solution)
        @warn("Preloaded solution not found, then it has been created\nPath: $path_solution")
    end

end

list_examples = Examples.list_examples

@testset "tests" begin

    for test in list_examples
        test_example(test.name, pem, test.d, test.N)
    end
end

