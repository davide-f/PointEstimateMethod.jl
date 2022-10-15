"""
    Distributions.moment(d::Normal, k)

Custom implementation of the moment function for Normal distributions
"""
function Distributions.moment(d::Normal, k)
    # central moment from https://en.wikipedia.org/wiki/Normal_distribution
    @assert d.μ == 0.0
    
    isodd(k) && return 0  # if odd return zero
    return (d.σ)^k * doublefactorial(k-1)  # return if even
end