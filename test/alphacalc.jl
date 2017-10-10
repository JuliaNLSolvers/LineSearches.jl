@testset "alpha-calculations" begin
    import NLSolversBase

    lsalphas =     [1.0,   0.5, 0.5, 0.49995, 0.5,  0.5]  # types
                    # Stat #HZ  wolfe   mt    bt2   bt3

    f(x) = vecdot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

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
end
