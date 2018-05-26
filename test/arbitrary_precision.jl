if Sys.WORD_SIZE == 32
    test_types = [Float32, Float64, BigFloat]
else
    test_types = [Float32, Float64, typeof(DoubleFloat64(1)), BigFloat]
end
for T in test_types
    @eval begin
        @testset "Arbitrary precision - initial step guess: $($T)" begin
            f(x) = vecdot(x, x)
            function g!(out, x)
                out[:] = 2x
            end

            x = $T[-1, -1]
            df = NLSolversBase.OnceDifferentiable(f,g!,x)

            phi_0, grtmp = NLSolversBase.value_gradient!(df, x)
            @test !isnan(phi_0)
            @test phi_0 isa $T
            @test !any(isnan, grtmp)
            @test grtmp isa Vector{$T}
            p = -grtmp
            dphi_0 = dot(p, grtmp)
            @test !isnan(dphi_0)
            @test dphi_0 isa $T

            function getstate()
                state = StateDummy($T(1),  x, similar(x), $T(NaN), p, $T(NaN))
            end
            # Test HagerZhang I0
            ls = HagerZhang{$T}()
            state = getstate()
            is = InitialHagerZhang{$T}(α0 = $T(NaN))
            is(ls, state, phi_0, dphi_0, df)
            @test !isnan(state.alpha)
            @test state.alpha isa $T
            @test ls.mayterminate[] == false

            # Test HagerZhang I12
            ls = HagerZhang{$T}()
            state = getstate()
            is = InitialHagerZhang{$T}(α0 = $T(1))
            is(ls, state, phi_0, dphi_0, df)
            @test !isnan(state.alpha)
            @test state.alpha isa $T
            @test ls.mayterminate[] == true

            # Test Static unscaled
            ls = HagerZhang{$T}()
            state = getstate()
            is = InitialStatic{$T}()
            is(ls, state, phi_0, dphi_0, df)
            @test state.alpha == is.alpha
            @test !isnan(state.alpha)
            @test state.alpha isa $T
            @test ls.mayterminate[] == false

            # Test Static scaled
            ls = HagerZhang{$T}()
            state = getstate()
            is = InitialStatic{$T}(alpha = $T(0.5), scaled = true)
            is(ls, state, phi_0, dphi_0, df)
            @test !isnan(state.alpha)
            @test !isnan(state.alpha)
            @test state.alpha isa $T
            @test ls.mayterminate[] == false

            # Test Previous
            ls = HagerZhang{$T}()
            state = getstate()
            alpha = state.alpha
            ls.mayterminate[] = true
            is = InitialPrevious{$T}()
            is(ls, state, phi_0, dphi_0, df)
            @test state.alpha == alpha
            @test !isnan(state.alpha)
            @test state.alpha isa $T
            @test ls.mayterminate[] == true

            # Test Previous NaN
            ls = HagerZhang{$T}()
            state = getstate()
            state.alpha = $T(NaN)
            ls.mayterminate[] = true
            is = InitialPrevious{$T}()
            is(ls, state, phi_0, dphi_0, df)
            @test state.alpha == is.alpha
            @test state.alpha isa $T
            @test ls.mayterminate[] == true

            # Test Quadratic NaN
            ls = HagerZhang{$T}()
            state = getstate()
            is = InitialQuadratic{$T}()
            is(ls, state, phi_0, dphi_0, df)
            @test state.alpha == is.α0
            @test !isnan(state.alpha)
            @test ls.mayterminate[] == false

            # Test Quadratic
            ls = HagerZhang{$T}()
            state = getstate()
            state.f_x_previous = 2*phi_0
            is = InitialQuadratic{$T}(snap2one=($T(0.9),$T(Inf)))
            is(ls, state, phi_0, dphi_0, df)
            @test !isnan(state.alpha)
            @test ls.mayterminate[] == false

            # Test Quadratic snap2one
            ls = HagerZhang{$T}()
            state = getstate()
            state.f_x_previous = 2*phi_0
            is = InitialQuadratic{$T}(snap2one=($T(0.75),$T(Inf)))
            is(ls, state, phi_0, dphi_0, df)
            @test !isnan(state.alpha)
            @test ls.mayterminate[] == false

            # Test ConstantChange NaN
            ls = HagerZhang{$T}()
            state = getstate()
            is = InitialConstantChange{$T}()
            is(ls, state, phi_0, dphi_0, df)
            @test state.alpha == is.α0
            @test !isnan(state.alpha)
            @test ls.mayterminate[] == false

            # Test ConstantChange
            ls = HagerZhang{$T}()
            state = getstate()
            state.dphi_0_previous = $T(0.1)*dphi_0
            is = InitialConstantChange{$T}()
            is(ls, state, phi_0, dphi_0, df)
            @test !isnan(state.alpha)
            @test ls.mayterminate[] == false

            # Test ConstantChange snap2one
            ls = HagerZhang{$T}()
            state = getstate()
            state.dphi_0_previous = $T(0.1)*dphi_0
            is = InitialConstantChange{$T}(snap2one=($T(0.25),$T(1)))
            is(ls, state, phi_0, dphi_0, df)
            @test !isnan(state.alpha)
            @test ls.mayterminate[] == false
        end
    end
end
