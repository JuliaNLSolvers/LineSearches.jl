# TODO: move the import into the let block when we stop supporting Julia v0.4?
import Optim

let
    println("\n#####################\nTest f_calls, g_calls")
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
    g!(x,out) = begin
        gcounter()
        prob.g!(x,out)
    end

    opts = Optim.Options(allow_f_increases=true)
    for ls in lsfunctions
        println("\nTesting $(string(ls))")
        fcounter(true); gcounter(true)

        res = Optim.optimize(f,g!, prob.h!, prob.initial_x,
                             Optim.BFGS(linesearch = ls), opts)
        @assert fcount == Optim.f_calls(res)
        @assert gcount == Optim.g_calls(res)
    end
end
