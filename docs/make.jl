using PointEstimateMethod
using Documenter

makedocs(
    modules = [PointEstimateMethod],
    doctest  = false,
    clean    = true,
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    sitename = "PointEstimateMethod.jl",
    authors  = "Davide Fioriti",
    pages   = [
        "Introduction" => [
            "index.md",
            "installation.md",
        ],
        "Examples" => [
            "examples/normal_distribution.md",
            "examples/general_distribution.md",
            "examples/multivariate_distribution.md",
        ],
        "API Reference" => "api_reference.md",
    ]
)

deploydocs(
    repo = "github.com/davide-f/PointEstimateMethod.jl.git",
    push_preview = true,
)
