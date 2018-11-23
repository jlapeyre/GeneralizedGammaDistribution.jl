using GeneralizedGammaDistribution
using Test

@testset "generalized gamma" begin
    @test GeneralizedGamma() == GeneralizedGamma(1.0, 1.0, 1.0)
end
