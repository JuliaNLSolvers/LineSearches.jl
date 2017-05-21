@testset "Optim usage" begin
    import Optim

    for ls in lsfunctions
        opts = Optim.Options(allow_f_increases = true)

        debug_printing && println("\nTesting $(string(ls))")
        counter = 0
        for (name, prob) in Optim.UnconstrainedProblems.examples
            if prob.isdifferentiable
                counter +=1
                if counter != 4
                    continue
                end

                res = Optim.optimize(prob.f, prob.g!, prob.h!,
                                     prob.initial_x, Optim.Newton(linesearch = ls), opts)
                debug_printing && println("$(name):\titerations = $(res.iterations)\tf_calls = $(Optim.f_calls(res))\tg_calls = $(Optim.g_calls(res))")
                @test Optim.minimum(res) <  prob.f(prob.solutions) + 1e-2
            end
        end
    end
end
