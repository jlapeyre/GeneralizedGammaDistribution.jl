module GenGammaDist
using Distributions

# Generalized Gamma distribution
# We have constructors for three different parameterizations.
# Some parameters can be negative, so this works as an inverse gamma distribution as well.

export GenGamma, gengamma_wiki, gengamma1, params, params1, params_wiki

# Use parameterization of R flexsurv implementation
immutable GenGamma
    μ
    σ
    Q
    gdist  # gamma distribution used to compute random samples
end

# construct with parameterization of R flexsurv implementation
function GenGamma(μ, σ, Q)
    gdist = Distributions.Gamma(1/Q^2,1)
    GenGamma(μ, σ, Q, gdist)
end

# Parameterization from wikipedia page
# PDF is: \frac{p/a^d}{\Gamma(d/p)} x^{d-1}e^{-(x/a)^p}
# The R page for gengamma flexsurv says that Q < 0 is allowed.
# The wikipedia page does not have abs(d) and abs(p), but this works, because the -i π cancels
# Must have d>0,p>0  or  d<0, p<0
function gengamma_wiki(a,d,p)
    μ = log(a) + (log(abs(d)) - log(abs(p)))/p
    σ = 1 / sqrt(p*d)
    Q = sqrt(p/d) * sign(p)
    GenGamma(μ,σ,Q)
end

# Another parameterization
# PDF is: \frac{p b^(d/p)}{\Gamma(d/p)} x^(d-1) e^{-b x^p}
# Must have d>0,p>0  or  d<0, p<0
function gengamma1(b,d,p)
    a = b^(-1/p)
    gengamma_wiki(a,d,p)
end



# Use the algorithm given on the R page in flexsurv
function Base.Random.rand(p::GenGamma)
    Qs = p.Q^2
#    gamma_deviate = rand(Distributions.Gamma(1/Qs,1)) # constructing all this is slow ?
    gamma_deviate = rand(p.gdist)  # only saves 10 or so percent time.
    w = log(Qs*gamma_deviate)/p.Q
    x = exp(p.μ + p.σ * w)
    return x
end

#  mean in wikipedia parameterization is
#  a * gamma((d+1)/p)/gamma(d/p)
function Base.mean(p::GenGamma)
    Qs = p.Q^2
    iQs = 1/Qs
    a = Qs^(p.σ / p.Q) * exp(p.μ)
    a * gamma(iQs + p.σ/p.Q)/gamma(iQs)
end


function params_wiki(dt::GenGamma)
    d = 1/(dt.σ * dt.Q)
    p = (dt.Q)/(dt.σ)
    a = abs(dt.Q)^(2/p)*exp(dt.μ)
    return (a,d,p)
end

function params1(dt::GenGamma)
    (a,d,p) = params_wiki(dt)
    b = a^(-p)
    return (b,d,p)
end

function params(d::GenGamma)
    return (d.μ, d.σ, d.Q)
end
end # module
