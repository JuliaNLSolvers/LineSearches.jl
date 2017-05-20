# This test checks the functionality cubic interpolation functionality
# in backtracking!
# The Himmelblau function converges, but require multiple iterations in
# some of the linesearches.

let
    import Optim

    name = "Himmelblau"
    prob = Optim.UnconstrainedProblems.examples[name]
    res = Optim.optimize(prob.f, prob.g!, prob.initial_x, Optim.BFGS(linesearch=LineSearches.BackTracking(order = 3)))
    @assert Optim.minimum(res) < prob.f(prob.solutions) + 1e-2
end
