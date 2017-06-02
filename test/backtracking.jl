# This test checks the functionality cubic interpolation functionality
# in BackTracking()
# The Himmelblau function converges, but require multiple iterations in
# some of the linesearches.

@testset "Check cubic backtracking" begin
    import Optim
    ls = LineSearches.BackTracking(order = 3)
    name = "Himmelblau"
    prob = Optim.UnconstrainedProblems.examples[name]
    res = Optim.optimize(prob.f, prob.g!, prob.initial_x, Optim.BFGS(linesearch=ls))
    @test Optim.minimum(res) < prob.f(prob.solutions) + 1e-2
end
