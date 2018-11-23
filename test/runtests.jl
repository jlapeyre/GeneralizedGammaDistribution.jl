using GenGammaDist
using Test

@testset "generalized gamma" begin
    GenGamma() == GenGamma(1.0, 1.0, 1.0)
end
