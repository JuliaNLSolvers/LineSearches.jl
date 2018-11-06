@testset "Initial step guess" begin
    f(x) = dot(x, x)
    function g!(out, x)
        out[:] = 2x
    end

    x = [-1., -1.]
    df = NLSolversBase.OnceDifferentiable(f,g!,x)

    phi_0, grtmp = NLSolversBase.value_gradient!(df, x)
    p = -grtmp
    dphi_0 = dot(p, grtmp)

    function getstate()
        state = StateDummy(1.0,  x, similar(x), NaN, p)
    end
    # Test HagerZhang I0
    ls = HagerZhang()
    state = getstate()
    is = InitialHagerZhang(α0 = NaN)
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == 0.005
    @test ls.mayterminate[] == false

    # Test HagerZhang I12
    ls = HagerZhang()
    state = getstate()
    is = InitialHagerZhang(α0 = 1.0)
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == 0.4999999999999999
    @test ls.mayterminate[] == true
    dfinf = NLSolversBase.NonDifferentiable(x->Inf, x)
    is = InitialHagerZhang(α0 = 1.0)
    is(ls, state, phi_0, dphi_0, dfinf)
    @test state.alpha == 0
    @test ls.mayterminate[] == true

    # Test Static unscaled
    ls = HagerZhang()
    state = getstate()
    is = InitialStatic()
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == is.alpha
    @test ls.mayterminate[] == false

    # Test Static scaled
    ls = HagerZhang()
    state = getstate()
    is = InitialStatic(alpha = 0.5, scaled = true)
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == 0.5 / norm(state.s)
    @test ls.mayterminate[] == false
    is = InitialStatic(alpha = 0.5, scaled = true)
    state.s .= (is.alpha / 100)
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha ==  1.0
    @test ls.mayterminate[] == false

    # Test Previous
    ls = HagerZhang()
    state = getstate()
    alpha = state.alpha
    ls.mayterminate[] = true
    is = InitialPrevious()
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == alpha
    @test ls.mayterminate[] == true

    # Test Previous NaN
    ls = HagerZhang()
    state = getstate()
    state.alpha = NaN
    ls.mayterminate[] = true
    is = InitialPrevious()
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == is.alpha
    @test ls.mayterminate[] == true # InitialPrevious should not touch this!

    # Test Quadratic NaN
    ls = HagerZhang()
    state = getstate()
    is = InitialQuadratic()
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == is.α0
    @test ls.mayterminate[] == false

    # Test Quadratic
    ls = HagerZhang()
    state = getstate()
    state.f_x_previous = 2*phi_0
    is = InitialQuadratic(snap2one=(0.9,Inf))
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == 1.0
    @test ls.mayterminate[] == false

    # Test Quadratic snap2one
    ls = HagerZhang()
    state = getstate()
    state.f_x_previous = 2*phi_0
    is = InitialQuadratic(snap2one=(0.75,Inf))
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == 1.0
    @test ls.mayterminate[] == false

    # Test ConstantChange NaN
    ls = HagerZhang()
    state = getstate()
    is = InitialConstantChange()
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == is.α0
    @test ls.mayterminate[] == false

    # Test ConstantChange
    ls = HagerZhang()
    state = getstate()
    is = InitialConstantChange()
    is.dϕ_0_previous[] = 0.1*dphi_0
    is(ls, state, phi_0, dphi_0, df)
    @test state.alpha == 0.25
    @test ls.mayterminate[] == false

    # Test ConstantChange snap2one
    ls = HagerZhang()
    state = getstate()
    is = InitialConstantChange(snap2one=(0.25,1.0))
    is.dϕ_0_previous[] = 0.1*dphi_0
    is(ls, state, phi_0, dphi_0, df)
    @test is.dϕ_0_previous[] == dphi_0
    @test state.alpha == 1.0
    @test ls.mayterminate[] == false
    dϕ_0_rand = rand()
    is(ls, state, phi_0, dϕ_0_rand, df)
    @test is.dϕ_0_previous[] == dϕ_0_rand
end
