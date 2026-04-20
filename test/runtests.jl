using Test
using Distributions
using YAML
using PointEstimateMethod

include("Examples.jl")

BASE_FOLDER = dirname(dirname(pathof(PointEstimateMethod)))

ATOL_TEST = 1e-4
RTOL_TEST = 1e-4
SEED = 0

function vector_approx_test(x_test::Vector, x_validation::Vector, rtol=RTOL_TEST, atol=ATOL_TEST)
    length(x_test) != length(x_validation) && return false
    return all(isapprox(xt, xv; rtol=rtol, atol=atol) for (xt, xv) in zip(x_test, x_validation))
end

"Function to test examples"
function test_example(example_name, testing_function, args...)

    # for reproducibility
    Distributions.Random.seed!(SEED)
    
    # calculate simulations
    calc_solution = testing_function(args...)

    path_solution = joinpath(BASE_FOLDER, "test", "testcases", string(testing_function), example_name * ".yml")

    I_calc = sortperm(calc_solution.x)
    
    if isfile(path_solution)
        # if the file exists run tests
        proven_solution = YAML.load_file(path_solution)
        
        @test vector_approx_test(calc_solution.x[I_calc], proven_solution["x"])
        @test vector_approx_test(calc_solution.p[I_calc], proven_solution["p"])
    else
        # otherwise create the tests
        mkpath(dirname(path_solution))

        dict_calc_solution = Dict(
            "x"=>calc_solution.x[I_calc],
            "p"=>calc_solution.p[I_calc],
        )

        YAML.write_file(path_solution, dict_calc_solution)
        @warn("Preloaded solution not found, then it has been created\nPath: $path_solution")
    end

end

"Function to test multivariate examples"
function test_multivariate_example(example_name, testing_function, args...)

    # for reproducibility
    Distributions.Random.seed!(SEED)

    # calculate simulations
    calc_solution = testing_function(args...)

    path_solution = joinpath(BASE_FOLDER, "test", "testcases", string(testing_function) * "_multivariate", example_name * ".yml")

    # Sort columns lexicographically for deterministic comparison
    I_calc = sortperm(eachcol(calc_solution.x), by=col -> collect(col))

    x_sorted = calc_solution.x[:, I_calc]
    p_sorted = calc_solution.p[I_calc]

    if isfile(path_solution)
        # if the file exists run tests
        proven_solution = YAML.load_file(path_solution)

        proven_x = Matrix(hcat([Float64.(row) for row in proven_solution["x"]]...)')
        proven_p = Float64.(proven_solution["p"])

        @test size(x_sorted) == size(proven_x)
        @test vector_approx_test(vec(x_sorted), vec(proven_x))
        @test vector_approx_test(p_sorted, proven_p)
    else
        # otherwise create the tests
        mkpath(dirname(path_solution))

        dict_calc_solution = Dict(
            "x" => [x_sorted[k, :] for k in 1:size(x_sorted, 1)],
            "p" => p_sorted,
        )

        YAML.write_file(path_solution, dict_calc_solution)
        @warn("Preloaded solution not found, then it has been created\nPath: $path_solution")
    end

end

list_examples = Examples.list_examples
list_multivariate_examples = Examples.list_multivariate_examples

@testset "tests" begin

    for test in list_examples
        test_example(test.name, pem, test.d, test.N)
    end

    for test in list_multivariate_examples
        test_multivariate_example(test.name, pem, test.d, test.N)
    end
end

