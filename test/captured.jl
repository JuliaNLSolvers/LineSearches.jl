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
    @testset "bisect! early Wolfe exit" begin
        # Directly test that bisect! returns early when it evaluates a
        # Wolfe-satisfying midpoint. Set up arrays matching bisect! preconditions:
        #   slopes[ia] < 0, values[ia] <= phi_lim
        #   slopes[ib] < 0, values[ib] > phi_lim
        alphas = [0.0, 1.0]
        values = [0.0, 5.0]
        slopes = [-2.0, -0.5]
        phi_0, dphi_0 = 0.0, -2.0
        phi_lim = phi_0 + 1e-6 * abs(phi_0)  # ≈ 0
        delta, sigma = 0.1, 0.9

        # ϕdϕ returns a point satisfying Wolfe at the midpoint d=0.5:
        #   Wolfe1: δ·dφ₀ ≥ (φ(d)-φ₀)/d  →  -0.2 ≥ -0.5/0.5 = -1.0  ✓
        #   Wolfe2: dφ(d) ≥ σ·dφ₀          →  -0.5 ≥ -1.8             ✓
        eval_count = Ref(0)
        function ϕdϕ_mock(d)
            eval_count[] += 1
            return (-0.5, -0.5)
        end

        ia, ib, wolfe = LineSearches.bisect!(
            ϕdϕ_mock, alphas, values, slopes, 1, 2,
            phi_lim, phi_0, dphi_0, delta, sigma,
        )
        @test wolfe
        @test eval_count[] == 1  # exited after just one evaluation
        @test alphas[ib] == 0.5  # the midpoint
        @test values[ib] == -0.5
        @test LineSearches.satisfies_wolfe(
            alphas[ib], values[ib], slopes[ib],
            phi_0, dphi_0, phi_lim, delta, sigma,
        )
    end

    @testset "bisect! without Wolfe continues normally" begin
        # When the midpoint does NOT satisfy Wolfe, bisect! should not
        # return early. Use a point with dφ ≥ 0 but high φ (normal exit).
        alphas = [0.0, 1.0]
        values = [0.0, 5.0]
        slopes = [-2.0, -0.5]
        phi_0, dphi_0 = 0.0, -2.0
        phi_lim = phi_0 + 1e-6 * abs(phi_0)
        delta, sigma = 0.1, 0.9

        # φ=0.5 > phi_lim and dφ=0.1 ≥ 0: triggers normal gphi≥0 exit.
        # Does NOT satisfy Wolfe because φ is above phi_lim and
        # sufficient decrease fails: δ·dφ₀ = -0.2 < (0.5-0)/0.5 = 1.0.
        ϕdϕ_high(d) = (0.5, 0.1)

        ia, ib, wolfe = LineSearches.bisect!(
            ϕdϕ_high, alphas, values, slopes, 1, 2,
            phi_lim, phi_0, dphi_0, delta, sigma,
        )
        @test !wolfe  # normal exit, not Wolfe
    end

    @testset "Wolfe during bracket expansion" begin
        # exp(α)-5α has minimum near ln(5)≈1.609. With c=0.1, the bracket
        # expansion loop (B3) multiplies c by ρ=5 repeatedly. At c=0.5
        # the approximate Wolfe conditions are satisfied, so the new check
        # in the bracketing loop should fire and return before entering
        # the main secant2! loop.
        ϕ_steep(α) = exp(α) - 5α
        dϕ_steep(α) = exp(α) - 5
        ϕdϕ_steep(α) = (ϕ_steep(α), dϕ_steep(α))

        ϕ0 = ϕ_steep(0.0)
        dϕ0 = dϕ_steep(0.0)

        # Verify preconditions: initial c=0.1 must NOT satisfy Wolfe
        # (so we actually enter the bracket expansion loop)
        phi_lim = ϕ0 + 1e-6 * abs(ϕ0)
        @test !LineSearches.satisfies_wolfe(0.1, ϕ_steep(0.1), dϕ_steep(0.1),
            ϕ0, dϕ0, phi_lim, 0.1, 0.9)

        # Verify that the expanded point c=0.5 DOES satisfy Wolfe
        @test LineSearches.satisfies_wolfe(0.5, ϕ_steep(0.5), dϕ_steep(0.5),
            ϕ0, dϕ0, phi_lim, 0.1, 0.9)

        # Run and verify: should terminate with α=0.5
        cache = LineSearchCache{Float64}()
        hz = HagerZhang(; cache)
        α, val = hz(ϕ_steep, ϕdϕ_steep, 0.1, ϕ0, dϕ0)
        @test α == 0.5
        @test val == ϕ_steep(0.5)
        # Only 2 cached evals: initial c=0.1, then expansion to c=0.5
        # (plus the implicit α=0 entry)
        @test length(cache.alphas) == 3
    end
end

@testset "PR#185" begin
    ϕ(α) = -exp(α)
    dϕ(α) = -exp(α)
    ϕdϕ(α) = ϕ(α), dϕ(α)

    α0 = 1.0
    ϕ0 = ϕ(0.0)
    dϕ0 = dϕ(0.0)

    _error_msg = @test_throws LineSearchException (HagerZhang(linesearchmax=14))(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
    @test _error_msg.value.alpha > α0 # Should be something gigantic like 1e34
end
