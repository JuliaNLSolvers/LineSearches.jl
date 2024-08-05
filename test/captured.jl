using LineSearches
using NLSolversBase
using Test

if !isdefined(@__MODULE__, :LineSearchTestCase)
    include("TestCases.jl")
    using .TestCases
end

@testset "Capturing data" begin
    cache = LineSearchCache{Float64}()
    lsalgs =  (HagerZhang(; cache), StrongWolfe(; cache), MoreThuente(; cache),
               BackTracking(; cache), BackTracking(; order=2, cache) )
    ϕ(x) = (x - π)^4
    dϕ(x) = 4*(x-π)^3
    ϕdϕ(x) = (ϕ(x), dϕ(x))
    for ls in lsalgs
        α, val = ls(ϕ, dϕ, ϕdϕ, 10.0, ϕ(0), dϕ(0))
        @test α < 10
        @test length(cache.alphas) == length(cache.values) && length(cache.alphas) > 1
    end
end

# From PR#174
@testset "PR#174" begin
    tc = LineSearchTestCase(
        [0.0, 1.0, 5.0, 3.541670844449739],
        [3003.592409634743, 2962.0378569864743, 2891.4462095232184, 3000.9760725116876],
        [-22332.321416890798, -20423.214551925797, 11718.185026267562, -22286.821227217057]
    )
    fdf = OnceDifferentiable(tc)
    hz = HagerZhang()
    α, val = hz(fdf.f, fdf.fdf, 1.0, fdf.fdf(0.0)...)
    @test_broken val <= minimum(tc)
end
