# GenGammaDist

Generalized gamma distribution, and inverse generalized gamma distribution.

Examples

Construct a distribution with same parameterization as the wikipedia page for
[Generalized gamma distribution](https://en.wikipedia.org/wiki/Generalized_gamma_distribution).
```julia
d = gengamma_wiki(a,d,p)
```

The parameters must satisfy either `d>0, p>0` or `d<0, p<0`.

Random sample
```julia
rand(d)
```

Compute the mean (which does not exist for `-1<d<0`)
```julia
mean(d)
```

Construct a distribution with same parameterization as the wikipedia page
except `b=a^(-p)`
```julia
d = gengamma1(b,d,p)
```

Construct a distribution with alternate parameterization
```julia
 GenGamma(μ, σ, Q)
```

The parameters are related via
`μ = log(a) + (log(abs(d)) - log(abs(p)))/p`,
`σ = 1/sqrt(p*d)`, `Q = sqrt(p/d)`, and
`d= 1/(σ*Q)`,  `p= Q/σ`,
`a= abs(Q)^(2*Q/σ)*exp(μ)`


Retrieve parameters

```julia
params(d)
params_wiki(d)
params1(d)
```

For `p=-1` and `d<0`, the distribution is the
[Inverse gamma distribution](https://en.wikipedia.org/wiki/Inverse-gamma_distribution),
with `α = -d`.
