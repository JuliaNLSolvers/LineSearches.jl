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

wolfec1(ls::BackTracking) = ls.c_1
wolfec1(ls::StrongWolfe) = ls.c_1
wolfec1(ls::MoreThuente) = ls.f_tol
wolfec1(ls::HagerZhang{T}) where T = zero(T)  # HZ uses Wolfe condition, but not necessarily based at zero

function amongtypes(x, types)
    t = typeof(x)
    return any(t <: ty for ty in types)
end

@testset "Captured cases" begin
    lsalgs =  (HagerZhang(), StrongWolfe(), MoreThuente(),
               BackTracking(), BackTracking(; order=2) )
    # From PR#174
    @testset "PR#174" begin
        tc = LineSearchTestCase(
            [0.0, 1.0, 5.0, 3.541670844449739],
            [3003.592409634743, 2962.0378569864743, 2891.4462095232184, 3000.9760725116876],
            [-22332.321416890798, -20423.214551925797, 11718.185026267562, -22286.821227217057]
        )
        fdf = OnceDifferentiable(tc)
        for ls in lsalgs
            α, val = ls(fdf.f, fdf.df, fdf.fdf, 1.0, fdf.fdf(0.0)...)
            @test 0 < α <= 5
            @test val < fdf.f(0) + wolfec1(ls) * α * fdf.df(0)   # Armijo (first Wolfe) condition
            if amongtypes(ls, (StrongWolfe, MoreThuente, HagerZhang))  # these types do not just backtrack
                @test val <= tc.values[2]
            end
        end
    end
    @testset "Issue#175" begin
        tc = LineSearchTestCase(
            [0.0, 0.2, 0.1, 0.055223623837026156],
            [3.042968312396456, 3.1174112871667603, -3.5035848233450224, 0.5244246783151265],
            [-832.4270136930788, -505.3362249257043, 674.9478303586366, 738.3388472427769]
        )
        fdf = OnceDifferentiable(tc)
        for ls in lsalgs
            α, val = ls(fdf.f, fdf.df, fdf.fdf, 0.2, fdf.fdf(0.0)...)
            @test 0 < α < 0.2
            @test val < fdf.f(0) + wolfec1(ls) * α * fdf.df(0)
        end
    end
end
