@testset "x_new" begin
    # On return, `x_new` must hold `x + α*s` for the returned α, so callers can
    # adopt that exact point instead of recomputing it.
    quadratic = (x = [-1.0, -1.0], f = x -> dot(x, x), g! = (out, x) -> out .= 2 .* x)
    himmelblau = let pr = OptimTestProblems.UnconstrainedProblems.examples["Himmelblau"]
        (x = copy(pr.initial_x), f = pr.f, g! = pr.g!)
    end

    @testset "$(typeof(ls)), $name" for ls in lstypes,
        (name, pr) in (("quadratic", quadratic), ("Himmelblau", himmelblau))

        (; x, f, g!) = pr
        df = NLSolversBase.OnceDifferentiable(f, g!, x)
        ϕ_0, g_0 = NLSolversBase.value_gradient!(df, x)
        s = -g_0
        x_new = similar(x)

        α, _ = ls(df, x, s, 1.0, x_new, ϕ_0, dot(s, g_0))

        @test x_new == muladd.(α, s, x)
    end
end
