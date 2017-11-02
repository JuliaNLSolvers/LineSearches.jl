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

    function getstate()
        lsr = LineSearchResults(eltype(x))
        push!(lsr, zero(phi0), phi0, dphi0)
        state = StateDummy(1.0,  x, similar(x), NaN, p, lsr, false)
    end

    # Test HagerZhang I0
    state = getstate()
    is = InitialHagerZhang()
    is(state, dphi0, df)
    @test state.alpha == 0.005
    @test state.mayterminate == true

    # Test HagerZhang I12
    state = getstate()
    is = InitialHagerZhang(α0 = 1.0)
    is(state, dphi0, df)
    @test state.alpha == 0.4999999999999999
    @test state.mayterminate == true

    # Test Static unscaled
    state = getstate()
    is = InitialStatic()
    is(state, dphi0, df)
    @test state.alpha == is.alpha
    @test state.mayterminate == true

    # Test Static scaled
    state = getstate()
    is = InitialStatic(alpha = 0.5, scaled = true)
    is(state, dphi0, df)
    @test state.alpha == 0.08838834764831843
    @test state.mayterminate == true

    # Test Previous
    state = getstate()
    state.mayterminate = false
    alpha = state.alpha
    is = InitialPrevious()
    is(state, dphi0, df)
    @test state.alpha == alpha
    @test state.mayterminate == false

    # Test Previous NaN
    state = getstate()
    state.alpha = NaN
    state.mayterminate = false
    is = InitialPrevious()
    is(state, dphi0, df)
    @test state.alpha == is.alpha
    @test state.mayterminate == true
end
