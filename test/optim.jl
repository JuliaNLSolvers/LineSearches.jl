# TODO: move the import into the let block when we stop supporting Julia v0.4?
import Optim

let
    for ls in lsfunctions
        println("\nTesting $(string(ls))")
        for (name, prob) in Optim.UnconstrainedProblems.examples
            if prob.isdifferentiable
                f_prob = prob.f
                res = Optim.optimize(f_prob, prob.initial_x, Optim.Newton(linesearch! = ls),
                                     Optim.OptimizationOptions(autodiff = true))
                println("$(name):\tf_calls = $(Optim.f_calls(res))\tg_calls = $(Optim.g_calls(res))")
                @assert norm(Optim.minimizer(res) - prob.solutions) < 1e-2
            end
        end
    end
end
