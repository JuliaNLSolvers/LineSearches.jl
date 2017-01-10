# TODO: move the import into the let block when we stop supporting Julia v0.4?
import Optim

let
    for ls in lsfunctions
        println("\nTesting $(string(ls))")
        for (name, prob) in Optim.UnconstrainedProblems.examples
            if prob.isdifferentiable
                res = Optim.optimize(prob.f, prob.initial_x, Optim.Newton(linesearch = ls),
                                     Optim.Options(autodiff = true))
                println("$(name):\titerations = $(res.iterations)\tf_calls = $(Optim.f_calls(res))\tg_calls = $(Optim.g_calls(res))")
                @assert Optim.minimum(res) <  prob.f(prob.solutions) + 1e-2
            end
        end
    end
end
