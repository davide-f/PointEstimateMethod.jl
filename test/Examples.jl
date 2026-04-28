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
        Example("Normal_truncated_3", truncated(Normal(1.0, 0.4), 0.0, +Inf), 3),
        Example("Normal_truncated_9", truncated(Normal(1.0, 0.4), 0.0, +Inf), 9),
    ]

    struct MultivariateExample
        name  # name of test
        d  # Vector of distributions
        N  # Number of estimate points (Integer or Vector)
    end

    list_multivariate_examples = [
        MultivariateExample("Normal2_N2", [Normal(), Normal()], 2),
        MultivariateExample("Normal2_N3", [Normal(), Normal()], 3),
        MultivariateExample(
            "Normal_TruncNormal_N2",
            [Normal(), truncated(Normal(1.0, 0.4), 0.0, +Inf)],
            2,
        ),
        MultivariateExample(
            "Normal_TruncNormal_mixed_N",
            [Normal(), truncated(Normal(1.0, 0.4), 0.0, +Inf)],
            [2, 3],
        ),
    ]

end