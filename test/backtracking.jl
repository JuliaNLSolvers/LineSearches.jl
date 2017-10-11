# This test checks the functionality cubic interpolation functionality
# in BackTracking()
# The Himmelblau function converges, but require multiple iterations in
# some of the linesearches.

@testset "Check cubic backtracking" begin
    pr = OptimTestProblems.Unconstrained.examples["Himmelblau"]
    x0 = copy(pr.initial_x)

    df = NLSolversBase.OnceDifferentiable(pr.f, pr.g!, x0)

    s = [42.0,18.0]
    lsr = LineSearchResults(eltype(x0))
    push!(lsr, 0.0, 26.0, -2088.0)

    mayterminate = false; alpha = 1.0; xtmp = zeros(x0)

    ls = BackTracking(order = 3)
    stepsize = ls(df, x0, s, xtmp, lsr, alpha, mayterminate)

    @test stepsize ≈ 0.020545340808876406
end
