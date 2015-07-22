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

Compute mean (which may not exist for negative parameters)
```julia
mean(d)
```

Construct a distribution with same parameterization as the wikipedia page
except `b=a^(-p)`
```julia
d = gengamma1(b,d,p)
```

For `p=-1` and `d<0`, the distribution is the
[Inverse gamma distribution](https://en.wikipedia.org/wiki/Inverse-gamma_distribution).

