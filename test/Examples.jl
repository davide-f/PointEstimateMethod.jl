module Examples

    using Distributions

    struct Example
        name  # name of test
        d  # Distributions
        N  # Number of estimate points
    end

    list_examples = [
        Example("Normal_2", Normal(), 2),
        Example("Normal_3", Normal(), 3),
    ]

end