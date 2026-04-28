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
        Example("MultiVariate_Normal_2", [Normal()], 2),
        Example("MultiVariate_2Normal_3", [Normal(), Normal()], 3),
        Example("MultiVariate_3Normal_2", [Normal(), Normal(), Normal()], 2),
        Example("MultiVariate_3Normal_3", [Normal(), Normal(), Normal()], 3),
    ]

end