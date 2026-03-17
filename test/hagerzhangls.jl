@testset "HagerZhangLS constructor validation" begin
    @test_throws ArgumentError HagerZhangLS(decrease = 0.0)
    @test_throws ArgumentError HagerZhangLS(decrease = -0.1)
    @test_throws ArgumentError HagerZhangLS(decrease = 0.5, curvature = 0.3)
    @test_throws ArgumentError HagerZhangLS(decrease = 0.6)
    @test_throws ArgumentError HagerZhangLS(curvature = 1.0)
    @test_throws ArgumentError HagerZhangLS(curvature = 1.5)
end

@testset "HagerZhangLS find_steplength" begin
    import LineSearches: find_steplength

    @testset "Simple quadratic: (α-2)²" begin
        φ(α) = ((α - 2)^2, 2(α - 2))
        φ0, dφ0 = φ(0.0)
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ, φ0, dφ0, 1.0)
        @test success
        @test α > 0
        @test fα < φ0
    end

    @testset "Quartic: (α-π)⁴" begin
        φ(α) = ((α - π)^4, 4(α - π)^3)
        φ0, dφ0 = φ(0.0)
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ, φ0, dφ0, 10.0)
        @test success
        @test fα < φ0
    end

    @testset "Rosenbrock 1D slice" begin
        φ(α) = begin
            v = (1 - α)^2 + 100*(α^2 - α)^2
            dv = -2*(1 - α) + 200*(α^2 - α)*(2α - 1)
            (v, dv)
        end
        φ0, dφ0 = φ(0.0)
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ, φ0, dφ0, 0.01)
        @test success
        @test fα < φ0
    end

    @testset "Non-finite initial step triggers backtracking" begin
        # φ blows up for α > 5
        φ(α) = α > 5 ? (Inf, Inf) : ((α - 2)^2, 2(α - 2))
        φ0, dφ0 = φ(0.0)
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ, φ0, dφ0, 100.0)
        @test success
        @test α ≤ 5
        @test isfinite(fα)
    end

    @testset "Non-finite φ0 returns failure" begin
        φ(α) = ((α - 2)^2, 2(α - 2))
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ, NaN, -4.0, 1.0)
        @test !success
        @test isnan(α)
    end

    @testset "Wolfe conditions satisfied at initial c" begin
        φ(α) = ((α - 2)^2, 2(α - 2))
        φ0, dφ0 = φ(0.0)
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ, φ0, dφ0, 2.0)
        @test success
        @test fα ≈ 0.0 atol=1e-10
    end

    @testset "Custom parameters" begin
        φ(α) = ((α - 2)^2, 2(α - 2))
        φ0, dφ0 = φ(0.0)
        hzl = HagerZhangLS(decrease=0.01, curvature=0.1, maxiter=100)
        α, fα, success = find_steplength(hzl, φ, φ0, dφ0, 1.0)
        @test success
        @test fα ≈ 0.0 atol=1e-6 # tight tolerances should get is near the optimum
    end

    @testset "Float32" begin
        φ(α) = ((α - 2)^2, 2(α - 2))
        φ0, dφ0 = φ(0.0f0)
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ, φ0, dφ0, 1.0f0)
        @test success
        @test α isa Float32
        @test fα isa Float32
    end

    @testset "Himmelblau via closures" begin
        pr = OptimTestProblems.UnconstrainedProblems.examples["Himmelblau"]
        x0 = copy(pr.initial_x)
        s = [42.0, 18.0]
        g = similar(x0)

        function φ_himmelblau(α)
            x = x0 .+ α .* s
            pr.g!(g, x)
            return (pr.f(x), dot(g, s))
        end

        φ0, dφ0 = φ_himmelblau(0.0)
        hzl = HagerZhangLS()
        α, fα, success = find_steplength(hzl, φ_himmelblau, φ0, dφ0, 1.0)
        @test success
        @test fα < φ0
    end
end
