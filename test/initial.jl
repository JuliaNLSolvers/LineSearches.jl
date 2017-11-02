@testset "Initial step guess" begin
    f(x) = vecdot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

    x = [-1., -1.]
    df = NLSolversBase.OnceDifferentiable(f,g!,x)

    phi0 = NLSolversBase.value_gradient!(df, x)
    grtmp = NLSolversBase.gradient(df)
    p = -grtmp
    dphi0 = dot(p, grtmp)

    lsr = LineSearchResults(eltype(x))
    push!(lsr, zero(phi0), phi0, dphi0)

    state = StateDummy(1.0,  x, similar(x), NaN, p, lsr, false)

    is = InitialHagerZhang()
    is(state, dphi0, df)
    @test state.alpha == 0.005
    @test state.mayterminate == true

    is = InitialHagerZhang(α0 = 1.0)
    is(state, dphi0, df)
    @test state.alpha == 0.49999999999994404
    @test state.mayterminate == true
end
