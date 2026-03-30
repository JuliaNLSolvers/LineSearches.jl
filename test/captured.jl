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
    @test val <= minimum(tc)
end

@testset "Early Wolfe termination" begin
    # Test that returned values satisfy Wolfe conditions and that
    # the early checks in bisect!/update!/bracket can reduce evaluations.

    # Helper: wrap ϕdϕ to count evaluations
    function counting_ϕdϕ(ϕdϕ_inner)
        count = Ref(0)
        function wrapper(α)
            count[] += 1
            ϕdϕ_inner(α)
        end
        return wrapper, count
    end

    # A function where the initial step overshoots, forcing bracket expansion
    # and bisection. The minimum is near α ≈ 1, but with initial c=10 the
    # algorithm must bracket and bisect, giving opportunities for early Wolfe.
    ϕ_quad(α) = (α - 1.0)^2
    dϕ_quad(α) = 2*(α - 1.0)
    ϕdϕ_quad(α) = (ϕ_quad(α), dϕ_quad(α))

    @testset "Wolfe satisfied on result" begin
        cache = LineSearchCache{Float64}()
        hz = HagerZhang(; cache)
        ϕdϕ_counted, count = counting_ϕdϕ(ϕdϕ_quad)
        α, val = hz(ϕ_quad, ϕdϕ_counted, 10.0, ϕ_quad(0.0), dϕ_quad(0.0))
        # The result must satisfy Wolfe conditions
        ϕ0 = ϕ_quad(0.0)
        dϕ0 = dϕ_quad(0.0)
        ϵ = hz.epsilon
        phi_lim = ϕ0 + ϵ * abs(ϕ0)
        @test LineSearches.satisfies_wolfe(α, val, dϕ_quad(α), ϕ0, dϕ0, phi_lim, hz.delta, hz.sigma)
        @test val == ϕ_quad(α)
    end

    @testset "Wolfe satisfied for steep exponential + quadratic" begin
        # A case that forces bracket expansion: steep initial descent then rise
        ϕ_steep(α) = exp(α) - 5α
        dϕ_steep(α) = exp(α) - 5
        ϕdϕ_steep(α) = (ϕ_steep(α), dϕ_steep(α))

        cache = LineSearchCache{Float64}()
        hz = HagerZhang(; cache)
        ϕ0 = ϕ_steep(0.0)
        dϕ0 = dϕ_steep(0.0)
        ϕdϕ_counted, count = counting_ϕdϕ(ϕdϕ_steep)
        α, val = hz(ϕ_steep, ϕdϕ_counted, 0.1, ϕ0, dϕ0)

        ϵ = hz.epsilon
        phi_lim = ϕ0 + ϵ * abs(ϕ0)
        @test LineSearches.satisfies_wolfe(α, val, dϕ_steep(α), ϕ0, dϕ0, phi_lim, hz.delta, hz.sigma)
        @test val ≈ ϕ_steep(α)
    end

    @testset "Result value matches cache minimum" begin
        # Using the test case from PR#174 which previously had @test_broken for this
        tc = LineSearchTestCase(
            [0.0, 1.0, 5.0, 3.541670844449739],
            [3003.592409634743, 2962.0378569864743, 2891.4462095232184, 3000.9760725116876],
            [-22332.321416890798, -20423.214551925797, 11718.185026267562, -22286.821227217057]
        )
        fdf = OnceDifferentiable(tc)
        cache = LineSearchCache{Float64}()
        hz = HagerZhang(; cache, check_flatness=true)
        α, val = hz(fdf.f, fdf.fdf, 1.0, fdf.fdf(0.0)...)
        @test minimum(cache.values) == val
    end
end

@testset "PR#185" begin
    ϕ(α) = 1/(α+1)
    dϕ(α) = -1/(α+1)^2
    ϕdϕ(α) = ϕ(α), dϕ(α)

    α0 = 1.0
    ϕ0 = ϕ(0.0)
    dϕ0 = dϕ(0.0)

    _error_msg = @test_throws LineSearchException (HagerZhang())(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
    @test _error_msg.value.alpha > α0 # Should be something gigantic like 1e34
end
