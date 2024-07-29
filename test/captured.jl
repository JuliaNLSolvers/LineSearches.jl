using LineSearches
using NLSolversBase
using Test

if !isdefined(@__MODULE__, :LineSearchTestCase)
    include("TestCases.jl")
    using .TestCases
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
