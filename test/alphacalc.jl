@testset "alpha-calculations" begin
    import NLSolversBase

    f(x) = vecdot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

    @testset "Quadratic function" begin
        lsalphas =     [1.0,   0.5, 0.5, 0.49995, 0.5,  0.5]  # types
                        # Stat #HZ  wolfe   mt    bt2   bt3
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
                println(string(linesearch!))
                alpha = linesearch!(df, x, p, xtmp, lsr, alpha, mayterminate)
                @test alpha == 1.0 # Is this what we want for non-descent directions?
            end
        end
    end

    @testset "Himmelblau" begin
        # This should be a bit more difficult, so hopefully it hits more of the algorithm steps
        function himmelblau(x::Vector)
            return (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
        end

        function himmelblau_gradient!(storage::Vector, x::Vector)
            storage[1] = 4.0 * x[1]^3 + 4.0 * x[1] * x[2] -
                44.0 * x[1] + 2.0 * x[1] + 2.0 * x[2]^2 - 14.0
            storage[2] = 2.0 * x[1]^2 + 2.0 * x[2] - 22.0 +
                4.0 * x[1] * x[2] + 4.0 * x[2]^3 - 28.0 * x[2]
        end

        f = himmelblau
        g! = himmelblau_gradient!

        x0 = [2.0, 2.0]
        df = NLSolversBase.OnceDifferentiable(himmelblau,himmelblau_gradient!,x0)

        s = [42.0,18.0]

        mayterminate = false; alpha = 1.0; xtmp = zeros(x0)

        lsalphas =  [1.0, # Stat
                     0.01375274240750926, # HZ
                     0.020646100006834013, # Wolfe
                     0.0175892025844326, # MT
                     0.020545340808876406, # BT(2)
                     0.010000000000000002] # BT(3)
        for (i, ls) in enumerate(lstypes)
            debug_printing && println("Testing $(string(ls))")
            lsr = LineSearchResults(eltype(x0))
            push!(lsr, 0.0, 26.0, -2088.0)

            df = NLSolversBase.OnceDifferentiable(f,g!,x0)
            stepsize = ls(df, x0, s, xtmp, lsr, alpha, mayterminate)

            @test stepsize == lsalphas[i]
        end
    end
end
