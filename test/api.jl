@testset "API" begin
    phi0 = 0.5
    dphi0 = 1.0
    alpha = 0.5
    lsr = LineSearchResults(typeof(phi0))

    push!(lsr, alpha, phi0, dphi0)
    @test lsr.alpha[1] == alpha
    @test lsr.value[1] == phi0
    @test lsr.slope[1] == dphi0
    @test lsr.nfailures == 0

    clear!(lsr)
    @test isempty(lsr.alpha)
    @test isempty(lsr.value)
    @test isempty(lsr.slope)
end
