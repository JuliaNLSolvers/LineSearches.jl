@testset "API" begin
    f(x) = vecdot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

    x = [-1., -1.]

    debug_printing && println("Testing alphatry")

    df = NLSolversBase.OnceDifferentiable(f,g!,x)


    phi0 = NLSolversBase.value_gradient!(df, x)
    grtmp = NLSolversBase.gradient(df)
    p = -grtmp
    dphi0 = dot(p, grtmp)

    lsr = LineSearchResults(eltype(x))
    push!(lsr, 0.5, phi0, dphi0)
    @test lsr.alpha[1] == 0.5
    @test lsr.value[1] == phi0
    @test lsr.slope[1] == dphi0

    alpha = alphainit(NaN, x,  grtmp, phi0)
    @test alpha == 0.005

    xtmp = zeros(x)
    alpha, mayterminate = alphatry(alpha, df, x, p, xtmp, lsr)
    @test alpha == 0.49999999999994404
    @test mayterminate == true

    clear!(lsr)
    @test isempty(lsr.alpha)
    @test isempty(lsr.value)
    @test isempty(lsr.slope)
end
