@testset "alpha-calculations" begin
    import NLSolversBase

    f(x) = dot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

    @testset "Quadratic function" begin
        lsalphas =     [1.0,   0.5, 0.5, 0.49995, 0.5,  0.5]  # types
                        # Stat  #HZ  wolfe   mt    bt3   bt2
        x = [-1., -1.]

        for (i, linesearch!) in enumerate(lstypes)
            debug_printing && println("Testing $(string(linesearch!))")
            df = NLSolversBase.OnceDifferentiable(f,g!,x)

            xtmp = copy(x)
            phi_0, grtmp = NLSolversBase.value_gradient!(df, x)
            p = -grtmp
            dphi_0 = dot(p, grtmp)

            alpha = 1.0
            if linesearch! == HagerZhang()
                linesearch!.mayterminate[] = false
            end

            alpha, ϕα = linesearch!(df, x, p, alpha, xtmp, phi_0, dphi_0)
            #xnew = x + alpha*p

            @test alpha ≈ lsalphas[i]
        end


        ### Test what happens at failures
        # x = [0, 0] is the minimum of f
        x = [0.0, 0.0]

        for (i, linesearch!) in enumerate(lstypes)
            debug_printing && println("Testing $(string(linesearch!))")
            df = NLSolversBase.OnceDifferentiable(f,g!,x)

            xtmp = copy(x)
            phi_0, grtmp = NLSolversBase.value_gradient!(df, x)
            p = -grtmp
            dphi_0 = dot(p, grtmp)

            alpha = 1.0
            if linesearch! == HagerZhang()
                linesearch!.mayterminate[] = false
            end

            if typeof(linesearch!) <: Union{HagerZhang, MoreThuente}
                @test_throws LineSearchException linesearch!(df, x, p, alpha,
                                                       xtmp, phi_0, dphi_0)
                try
                    linesearch!(df, x, p, alpha,
                                xtmp, phi_0, dphi_0)
                catch ex
                    @test ex.message == "Search direction is not a direction of descent."
                end
            else
                α, ϕα = linesearch!(df, x, p, alpha, xtmp, phi_0, dphi_0)
                @test ϕα == 0.0
                @test alpha == α # Is this what we want for zero-slope directions?
            end
        end
    end

    @testset "Himmelblau" begin
        # This should be a bit more difficult, so hopefully it hits more of the algorithm steps
        pr = OptimTestProblems.UnconstrainedProblems.examples["Himmelblau"]
        x0 = copy(pr.initial_x)


        s = [42.0,18.0]

        mayterminate = Ref{Bool}(false)
        alpha = 1.0; xtmp = zero(x0)

        lsalphas =  [1.0,                  # Static
                     0.01375274240750926,  # HZ
                     0.020646100006834013, # Wolfe
                     0.0175892025844326,   # MT
                     0.020545340808876406, # BT(3)
                     0.010000000000000002] # BT(2)
        for (i, ls) in enumerate(lstypes)
            debug_printing && println("Testing $(string(ls))")

            phi_0 = 26.0
            dphi_0 = -2088.0

            df = NLSolversBase.OnceDifferentiable(pr.f, pr.g!, x0)

            stepsize, ϕstepsize = ls(df, x0, s, alpha, xtmp, phi_0, dphi_0)

            @test stepsize ≈ lsalphas[i]
        end
    end
end
