@testset "Counters" begin
    import Optim

    debug_printing && println("\n#####################\nTest f_calls, g_calls")
    prob = Optim.UnconstrainedProblems.examples["Rosenbrock"]

    let
        global fcount = 0
        global fcounter
        function fcounter(reset::Bool = false)
            if reset
                fcount = 0
            else
                fcount += 1
            end
            fcount
        end
        global gcount = 0
        global gcounter
        function gcounter(reset::Bool = false)
            if reset
                gcount = 0
            else
                gcount += 1
            end
            gcount
        end
    end

    f(x) = begin
        fcounter()
        prob.f(x)
    end
    g!(out, x) = begin
        gcounter()
        prob.g!(out, x)
    end

    for ls in lstypes
        debug_printing && println("\nTesting $(string(ls))")
        fcounter(true); gcounter(true)

        res = Optim.optimize(f,g!, prob.initial_x,
                             Optim.BFGS(linesearch = ls))
        # @test fcount == Optim.f_calls(res)
        # @test gcount == Optim.g_calls(res)
        @test fcount == Optim.f_calls(res) + 1
        @test gcount == Optim.g_calls(res) + 3
    end
end
