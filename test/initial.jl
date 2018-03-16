@testset "Initial step guess" begin
    f(x) = vecdot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

    x = [-1., -1.]
    df = NLSolversBase.OnceDifferentiable(f,g!,x)

    phi_0 = NLSolversBase.value_gradient!(df, x)
    grtmp = NLSolversBase.gradient(df)
    p = -grtmp
    dphi_0 = dot(p, grtmp)

    function getstate()
        state = StateDummy(1.0,  x, similar(x), NaN, p, Ref{Bool}(false), NaN)
    end

    # Test HagerZhang I0
    state = getstate()
    is = InitialHagerZhang(α0 = NaN)
    is(state, phi_0, dphi_0, df)
    @test state.alpha == 0.005
    @test state.mayterminate[] == false

    # Test HagerZhang I12
    state = getstate()
    is = InitialHagerZhang(α0 = 1.0)
    is(state, phi_0, dphi_0, df)
    @test state.alpha == 0.4999999999999999
    @test state.mayterminate[] == true

    # Test Static unscaled
    state = getstate()
    is = InitialStatic()
    is(state, phi_0, dphi_0, df)
    @test state.alpha == is.alpha
    @test state.mayterminate[] == false

    # Test Static scaled
    state = getstate()
    is = InitialStatic(alpha = 0.5, scaled = true)
    is(state, phi_0, dphi_0, df)
    @test state.alpha == 0.08838834764831843
    @test state.mayterminate[] == false

    # Test Previous
    state = getstate()
    alpha = state.alpha
    state.mayterminate[] = true
    is = InitialPrevious()
    is(state, phi_0, dphi_0, df)
    @test state.alpha == alpha
    @test state.mayterminate[] == true

    # Test Previous NaN
    state = getstate()
    state.alpha = NaN
    state.mayterminate[] = true
    is = InitialPrevious()
    is(state, phi_0, dphi_0, df)
    @test state.alpha == is.alpha
    @test state.mayterminate[] == true # InitialPrevious should not touch this!

    # Test Quadratic NaN
    state = getstate()
    is = InitialQuadratic()
    is(state, phi_0, dphi_0, df)
    @test state.alpha == is.α0
    @test state.mayterminate[] == false

    # Test Quadratic
    state = getstate()
    state.f_x_previous = 2*phi_0
    is = InitialQuadratic(snap2one=(0.9,Inf))
    is(state, phi_0, dphi_0, df)
    @test state.alpha == 0.8200000000000001
    @test state.mayterminate[] == false

    # Test Quadratic snap2one
    state = getstate()
    state.f_x_previous = 2*phi_0
    is = InitialQuadratic(snap2one=(0.75,Inf))
    is(state, phi_0, dphi_0, df)
    @test state.alpha == 1.0
    @test state.mayterminate[] == false

    # Test ConstantChange NaN
    state = getstate()
    is = InitialConstantChange()
    is(state, phi_0, dphi_0, df)
    @test state.alpha == is.α0
    @test state.mayterminate[] == false

    # Test ConstantChange
    state = getstate()
    state.dphi_0_previous = 0.1*dphi_0
    is = InitialConstantChange()
    is(state, phi_0, dphi_0, df)
    @test state.alpha == 0.25
    @test state.mayterminate[] == false

    # Test ConstantChange snap2one
    state = getstate()
    state.dphi_0_previous = 0.1*dphi_0
    is = InitialConstantChange(snap2one=(0.25,1.0))
    is(state, phi_0, dphi_0, df)
    @test state.alpha == 1.0
    @test state.mayterminate[] == false
end
