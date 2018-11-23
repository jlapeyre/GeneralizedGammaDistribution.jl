"""
    module GeneralizedGammaDistribution

This module provides the generalized gamma distribution.

Constructors: GeneralizedGamma

Parameterizations:

Three parameterizations are supported. They are represented by types `RflexsurvParams`, `StacyParams`, and `WikipediaParams`.
`RflexsurvParams` is the default parameterization.

These types may be given as the first argument to the constructor `GeneralizedGamma`, and to `params`.
"""
module GeneralizedGammaDistribution
import Distributions

export GeneralizedGamma
export RflexsurvParams, StacyParams, WikipediaParams

# Use parameterization of R flexsurv implementation
"""
    GeneralizedGamma(μ=1.0, σ=1.0, Q=1.0)

The *Generalized gamma distribution* has probability density
function (in the Wikipedia parameterization)

```math
f(x; a, d, p) = \\frac{(p/a^d) x^{d-1} e^{-(x/a)^p}}{\\Gamma(d/p)}
\\quad x > 0
```

```julia
params(RflexsurvParams, d)  # Get the parameters, i.e. (μ, σ, Q)
params(d)                   # Get the parameters, i.e. (μ, σ, Q)
params(WikipediaParams, d)  # Get the parameters, (a, d, p)
params(StacyParams, d)      # Get the parameters, (b, d, p)
```

External links

* [Generalized gamma distribution on Wikipedia](https://en.wikipedia.org/wiki/Generalized_gamma_distribution)
"""
struct GeneralizedGamma{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    μ::T
    σ::T
    Q::T
    gamma_dist::Distributions.Gamma  # gamma distribution used to compute random samples
end

"""
    RflexsurvParams

Parameterization of the generalized gamma distribution used in the R flexsurve implementation.
This is the default parameterization
"""
struct RflexsurvParams # <: AbstractParameterization
end

"""
    WikipediaParams

Parameterization of the generalized gamma distribution found in the Wikipedia article
"""
struct WikipediaParams # <: AbstractParameterization
end

"""
    StacyParams

The Stacy Parameterization of the generalized gamma distribution. This parameterization was
used in the older dgengamma.orig implementation.
"""
struct StacyParams # <: AbstractParameterization
end

"""
    GeneralizedGamma(::Type{RflexsurvParams}, μ, σ, Q)
    GeneralizedGamma(μ, σ, Q)

The `generalized gamma distribution` parameterized as in the R `flexsurv`
implementation[^1].

[^1]: Prentice, R. L. (1974). A log gamma model and its maximum likelihood estimation. Biometrika 61(3):539-544
"""
function GeneralizedGamma(::Type{RflexsurvParams}, μ=1.0, σ=1.0, Q=1.0)
    gamma_dist = Distributions.Gamma(1 / Q^2, 1)
    return GeneralizedGamma(μ, σ, Q, gamma_dist)
end
GeneralizedGamma(μ::Real=1.0, σ::Real=1.0, Q::Real=1.0) = GeneralizedGamma(RflexsurvParams, μ, σ, Q)

# Parameterization from wikipedia page
# PDF is: \frac{p/a^d}{\Gamma(d/p)} x^{d-1}e^{-(x/a)^p}
# The R page for gengamma flexsurv says that Q < 0 is allowed.
# The wikipedia page does not have abs(d) and abs(p), but this works, because the -i π cancels
# Must have d > 0, p > 0  or  d < 0, p < 0
"""
    GeneralizedGamma(::Type{WikipediaParams}, a=1.0, d=1.0, p=1.0)

The `generalized gamma distribution` parameterized as on the Wikipedia page.
Requires either `d > 0, p > 0` or `d < 0, p < 0`.
"""
function GeneralizedGamma(::Type{WikipediaParams}, a=1.0, d=1.0, p=1.0)
    μ = log(a) + (log(abs(d)) - log(abs(p))) / p
    σ = 1 / sqrt(p * d)
    Q = sqrt(p / d) * sign(p)
    return GeneralizedGamma(μ, σ, Q)
end

"""
    GeneralizedGamma(::Type{StacyParams}, b=1.0, d=1.0, p=1.0)

The `generalized gamma distribution` with the Stacy parameterization. This follows
the older R convention `dgengamma.orig` [^2].

[^2]: Stacy, E. W. (1962). A generalization of the gamma distribution. Annals of Mathematical Statistics 33:1187-92.
"""
function GeneralizedGamma(::Type{StacyParams}, b=1.0, d=1.0, p=1.0)
    a = b^(-1 / p)
    return GeneralizedGamma(WikipediaParams, a, d, p)
end

@inline Distributions.partype(d::GeneralizedGamma{T}) where {T<:Real} = T

# Use the algorithm given in the documentation pages for the R package flexsurv.
function Distributions.rand(p::GeneralizedGamma)
    Qs = p.Q^2
    gamma_deviate = rand(p.gamma_dist)  # only saves 10 or so percent time.
    w = log(Qs * gamma_deviate) / p.Q
    return exp(p.μ + p.σ * w)
end

#  mean in wikipedia parameterization is
#  a * gamma((d+1)/p)/gamma(d/p)
function Distributions.mean(p::GeneralizedGamma)
    Qs = p.Q^2
    iQs = 1 / Qs
    a = Qs^(p.σ / p.Q) * exp(p.μ)
    return a * Distributions.gamma(iQs + p.σ / p.Q) / Distributions.gamma(iQs)
end

function Distributions.pdf(p::GeneralizedGamma, x::Real)
    (d, p, a) = Distributions.params(WikipediaParams, p)
    return (p / a^d)/ Distributions.gamma(d / p) * x^(d - 1) * exp(-(x / a)^p)
end

function Distributions.mode(p::GeneralizedGamma{T}) where T
    (d, p, a) = Distributions.params(WikipediaParams, p)
    return d <= 1 ? zero(T) : a * ((d - 1) / p)^(1/p)
end

Distributions.support(::GeneralizedGamma{T}) where {T} = Distributions.RealInterval(zero(T), convert(T, Inf))

function Distributions.var(dt::GeneralizedGamma)
    (d, p, a) = Distributions.params(WikipediaParams, dt)
    d2 = Distributions.gamma(d / p)
    return a^2 * ((Distributions.gamma((d + 2) / p) / d2) -
                  (Distributions.gamma((d + 1) / p) / d2)^2)
end

# How do we include the lower incomplete gamma function ? How is it done in Distributions ?
# It is buried somewhere in the layers of code.
# function Distributions.cdf(p::GeneralizedGamma, x::Real)
#     (d, p, a) = Distributions.params(WikipediaParams, p)
#     num = lower_incomplete_gamma(d/p, (x/a)^p)
#     return num /  Distributions.gamma(d / p)
# end

"""
    params(::Type{WikipediaParams}, dt::GeneralizedGamma)

Return parameters of the `generalized Gamma distribution` following
the convention on the Wikipedia page.
"""
function Distributions.params(::Type{WikipediaParams}, dt::GeneralizedGamma)
    d = 1 / (dt.σ * dt.Q)
    p = (dt.Q) / (dt.σ)
    a = abs(dt.Q)^(2 / p) * exp(dt.μ)
    return (a, d, p)
end

"""
   params(::Type{StacyParams}, dt::GeneralizedGamma)

Return parameters of the `generalized Gamma distribution` following the parameterization
of the older R implementation `dgengamma.orig` [^2].

[^2]: Stacy, E. W. (1962). A generalization of the gamma distribution. Annals of Mathematical Statistics 33:1187-92.
"""
function Distributions.params(::Type{StacyParams}, dt::GeneralizedGamma)
    (a, d, p) = Distributions.params(WikipediaParams, dt)
    b = a^(-p)
    return (b, d, p)
end

"""
   params(::Type{RflexsurvParams}, dt::GeneralizedGamma)

Return parameters of the `generalized Gamma distribution` following
the convention of the R `flexsurv` function.
"""
function Distributions.params(::Type{RflexsurvParams}, d::GeneralizedGamma)
    return (d.μ, d.σ, d.Q)
end

Distributions.params(d::GeneralizedGamma) = Distributions.params(RflexsurvParams, d)

end # module
