import Optim

# This test checks the functionality cubic interpolation functionality
# in backtracking!
# The Himmelblau function converges, but require multiple iterations in
# some of the linesearches.

let
    name = "Himmelblau"
    prob = Optim.UnconstrainedProblems.examples[name]
    res = Optim.optimize(prob.f, prob.initial_x, Optim.BFGS(linesearch=LineSearches.bt3!),
                         Optim.Options(autodiff = true))

    @assert norm(Optim.minimizer(res) - prob.solutions) < 1e-2
end
