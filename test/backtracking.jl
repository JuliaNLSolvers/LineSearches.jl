# This test checks the functionality cubic interpolation functionality
# in BackTracking()
# The Himmelblau function converges, but require multiple iterations in
# some of the linesearches.

@testset "Check cubic backtracking" begin

    function himmelblau(x::Vector)
        return (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
    end

    function himmelblau_gradient!(storage::Vector, x::Vector)
        storage[1] = 4.0 * x[1]^3 + 4.0 * x[1] * x[2] -
            44.0 * x[1] + 2.0 * x[1] + 2.0 * x[2]^2 - 14.0
        storage[2] = 2.0 * x[1]^2 + 2.0 * x[2] - 22.0 +
            4.0 * x[1] * x[2] + 4.0 * x[2]^3 - 28.0 * x[2]
    end
    x0 = [2.0, 2.0]

    df = NLSolversBase.OnceDifferentiable(himmelblau,himmelblau_gradient!,x0)

    s = [42.0,18.0]
    lsr = LineSearchResults(eltype(x0))
    push!(lsr, 0.0, 26.0, -2088.0)

    mayterminate = false; alpha = 1.0; xtmp = zeros(x0)

    ls = BackTracking(order = 3)
    stepsize = ls(df, x0, s, xtmp, lsr, alpha, mayterminate)

    @test stepsize ≈ 0.020545340808876406
end
