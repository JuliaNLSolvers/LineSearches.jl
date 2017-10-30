@testset "alpha-calculations" begin
    import NLSolversBase

    f(x) = vecdot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

    @testset "Quadratic function" begin
        lsalphas =     [1.0,   0.35355339059327373, 0.5, 0.5, 0.49995, 0.5,  0.5]  # types
                        # Stat # Stat(scaled)       #HZ  wolfe   mt    bt2   bt3
        x = [-1., -1.]

        for (i, linesearch!) in enumerate(lstypes)
            debug_printing && println("Testing $(string(linesearch!))")
            df = NLSolversBase.OnceDifferentiable(f,g!,x)

            xtmp = copy(x)
            phi0 = NLSolversBase.value_gradient!(df, x)
            grtmp = NLSolversBase.gradient(df)
            p = -grtmp
            dphi0 = dot(p, grtmp)

            lsr = LineSearchResults(eltype(x))
            push!(lsr, 0.0, phi0, dphi0)

            alpha = 1.0
            mayterminate = false

            alpha = linesearch!(df, x, p, xtmp, lsr, alpha, mayterminate)
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
            phi0 = NLSolversBase.value_gradient!(df, x)
            grtmp = NLSolversBase.gradient(df)
            p = -grtmp
            dphi0 = dot(p, grtmp)

            lsr = LineSearchResults(eltype(x))
            push!(lsr, 0.0, phi0, dphi0)

            alpha = 1.0
            mayterminate = false

            if linesearch! == HagerZhang() || linesearch! == MoreThuente()
                @test_throws ErrorException alpha = linesearch!(df, x, p, xtmp, lsr, alpha, mayterminate)
            else
                alpha = linesearch!(df, x, p, xtmp, lsr, alpha, mayterminate)
                @test alpha == 1.0 # Is this what we want for non-descent directions?
            end
        end
    end

    @testset "Himmelblau" begin
        # This should be a bit more difficult, so hopefully it hits more of the algorithm steps
        pr = OptimTestProblems.UnconstrainedProblems.examples["Himmelblau"]
        x0 = copy(pr.initial_x)


        s = [42.0,18.0]

        mayterminate = false; alpha = 1.0; xtmp = zeros(x0)

        lsalphas =  [1.0,                  # Static(scaled)
                     0.021884405476620426, # Static(scaled)
                     0.01375274240750926,  # HZ
                     0.020646100006834013, # Wolfe
                     0.0175892025844326,   # MT
                     0.020545340808876406, # BT(2)
                     0.010000000000000002] # BT(3)
        for (i, ls) in enumerate(lstypes)
            debug_printing && println("Testing $(string(ls))")
            lsr = LineSearchResults(eltype(x0))
            push!(lsr, 0.0, 26.0, -2088.0)
            df = NLSolversBase.OnceDifferentiable(pr.f, pr.g!, x0)

            stepsize = ls(df, x0, s, xtmp, lsr, alpha, mayterminate)

            @test stepsize ≈ lsalphas[i]
        end
    end
end
