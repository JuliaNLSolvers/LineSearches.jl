# Provide functionality for initial guesses to the step length

"""
Provide static initial step length.

Keyword `alpha` corresponds to static step length, default is 1.0.
If keyword `scaled = true`, then the initial step length
is scaled with the `l_2` norm of the step direction.
"""
@with_kw struct InitialStatic{T}
    alpha::T = 1.0
    scaled::Bool = false # Scales step. alpha ← min(alpha,||s||_2) / ||s||_2
end

function (is::InitialStatic{T})(ls, state, phi_0, dphi_0, df) where T
    PT = promote_type(T, real(eltype(state.s)))
    if is.scaled == true && (ns = real(norm(state.s))) > convert(PT, 0)
        # TODO: Type instability if there's a type mismatch between is.alpha and ns?
        state.alpha = convert(PT, min(is.alpha, ns)) / ns
    else
        state.alpha = convert(PT, is.alpha)
    end
end

"""
Use previous step length as initial guess,
within the bounds [alphamin, alphamax]


If state.alpha is NaN, then return fallback value is.alpha
"""
@with_kw struct InitialPrevious{T}
    alpha::T = 1.0
    alphamin::T = 0.0
    alphamax::T = Inf
end

function (is::InitialPrevious)(ls, state, phi_0, dphi_0, df)
    if isnan(state.alpha)
        state.alpha = is.alpha
    end
    state.alpha = min(is.alphamax, state.alpha)
    state.alpha = max(is.alphamin, state.alpha)
end



"""
Quadratic interpolation for initial step length guess.

This is meant for methods that do not produce well-scaled search directions,
such as Gradient Descent and (variations of) Conjugate Gradient methods.
See the discussion around Nocedal and Wright, 2nd ed, (3.60).

This procedure have several arguments, with the following defaults.
- `α0       = 1.0`.         The initial step size at the first iteration.
- `αmin     = 1e-12`.       The minimum initial step size. (Default arbitrary).
- `αmax     = 1.0`.         The maximum initial step size.
- `ρ        = 0.25`.        Maximum decrease from previous iteration, `αinit ≥ α_{k-1}`. (Default arbitrary).
- `snap2one = (0.75, Inf)`. Set all values within this (closed) interval to 1.0. (Default arbitrary).

If αmax ≠ 1.0, then you should consider to ensure that snap2one[2] < αmax.
"""
@with_kw struct InitialQuadratic{T}
    αmin::T = 1e-12 # Minimum initial step size (value somewhat arbitrary)
    αmax::T = 1.0   # Maximum initial step size (advised by Nocedal+Wright)
    α0::T   = 1.0   # Fallback at first iteration
    ρ::T    = 0.25  # maximum decrease from previous step (value somewhat arbitrary)
    snap2one::Tuple{T,T} = (0.75, Inf) # Set everything in this (closed) interval to 1.0
end

function (is::InitialQuadratic{T})(ls, state, phi_0, dphi_0, df) where T
    if !isfinite(state.f_x_previous) || isapprox(dphi_0, convert(T, 0), atol=eps(T)) # Need to add a tolerance
        # If we're at the first iteration
        αguess = is.α0
    else
        αguess = 2 * (NLSolversBase.value(df) - state.f_x_previous) / dphi_0
        αguess = NaNMath.max(is.αmin, state.alpha*is.ρ, αguess)
        αguess = NaNMath.min(is.αmax, αguess)
        # if αguess ≈ 1, then make it 1 (Newton-type behaviour)
        if is.snap2one[1] ≤ αguess ≤ is.snap2one[2]
            αguess = one(state.alpha)
        end
    end
    state.alpha = αguess
end

"""
Constant first-order change approximation to determine initial step length.

** This requires that the optimization algorithm stores dphi0 from the previous iteration **
(dphi0_previous = real(dot(∇f_{k-1}, s_{k-1})), where s is the step direction.

This is meant for methods that do not produce well-scaled search directions,
such as Gradient Descent and (variations of) Conjugate Gradient methods.
See the discussion in Nocedal and Wright, 2nd ed, p. 59 on "Initial Step Length"

This procedure have several arguments, with the following defaults.
- `α0       = 1.0`.         The initial step size at the first iteration.
- `αmin     = 1e-12`.       The minimum initial step size. (Default arbitrary).
- `αmax     = 1.0`.         The maximum initial step size.
- `ρ        = 0.25`.        Maximum decrease from previous iteration, `αinit ≥ α_{k-1}`. (Default arbitrary).
- `snap2one = (0.75, Inf)`. Set all values within this (closed) interval to 1.0. (Default arbitrary).

If αmax ≠ 1.0, then you should consider to ensure that snap2one[2] < αmax.
"""
struct InitialConstantChange{T}
    αmin::T # Minimum initial step size (value somewhat arbitrary)
    αmax::T # Maximum initial step size
    α0::T   # Fallback at first iteration
    ρ::T    # maximum decrease from previous step
    snap2one::Tuple{T,T} # Set everything in this (closed) interval to 1.0
    dϕ_0_previous::Base.RefValue{T}
end

# Have to make this constructor without with_kw because Ref(NaN) has to adapt to T
function InitialConstantChange{T}(; αmin = 1e-12,
                        αmax = 1.0,
                        α0   = 1.0,
                        ρ    = 0.25,
                        snap2one = (0.75, Inf)) where T
    αmin, αmax, α0, ρ = convert.(T, (αmin, αmax, α0, ρ))
    snap2one = convert.(T, snap2one)
    InitialConstantChange(αmin, αmax, α0, ρ, snap2one, Ref{T}(convert(T, NaN)))
end

# Have to make this constructor without with_kw because Ref(NaN) has to adapt to T
function InitialConstantChange(; αmin = 1e-12,
                        αmax = 1.0,
                        α0   = 1.0,
                        ρ    = 0.25,
                        snap2one = (0.75, Inf))
    T = promote_type(typeof.((αmin, αmax, α0, ρ))...)
    InitialConstantChange(αmin, αmax, α0, ρ, snap2one, Ref{T}(convert(T, NaN)))
end

function (is::InitialConstantChange{T})(ls, state, phi_0, dphi_0, df) where T
    if !isfinite(is.dϕ_0_previous[]) || !isfinite(state.alpha) ||
        isapprox(dphi_0, convert(T, 0), atol=eps(T))
        # If we're at the first iteration
        αguess = is.α0
    else
        # state.alpha is the previously used step length
        αguess = state.alpha * is.dϕ_0_previous[] / dphi_0
        αguess = NaNMath.max(is.αmin, state.alpha*is.ρ, αguess)
        αguess = NaNMath.min(is.αmax, αguess)
        # if αguess ≈ 1, then make it 1 (Newton-type behaviour)
        if is.snap2one[1] ≤ αguess ≤ is.snap2one[2]
            αguess = one(state.alpha)
        end
    end
    is.dϕ_0_previous[] = dphi_0
    state.alpha = αguess
end


"""
Initial step size algorithm from
  W. W. Hager and H. Zhang (2006) Algorithm 851: CG_DESCENT, a
    conjugate gradient method with guaranteed descent. ACM
    Transactions on Mathematical Software 32: 113–137.

If α0 is NaN, then procedure I0 is called at the first iteration,
otherwise, we select according to procedure I1-2, with starting value α0.
"""
@with_kw struct InitialHagerZhang{T}
    ψ0::T          = 0.01
    ψ1::T          = 0.2
    ψ2::T          = 2.0
    ψ3::T          = 0.1
    αmax::T        = Inf
    α0::T          = 1.0 # Initial alpha guess. NaN => algorithm calculates
    quadstep::Bool = true
    verbose::Bool  = false
end

function (is::InitialHagerZhang)(ls::Tls, state, phi_0, dphi_0, df) where Tls
    if isnan(state.f_x_previous) && isnan(is.α0)
        # If we're at the first iteration (f_x_previous is NaN)
        # and the user has not provided an initial step size (is.α0 is NaN),
        # then we
        # pick the initial step size according to HZ #I0
        state.alpha = _hzI0(state.x, NLSolversBase.gradient(df),
                            NLSolversBase.value(df),
                            is.αmax,
                            convert(eltype(state.x), is.ψ0)) # Hack to deal with type instability between is{T} and state.x
        if Tls <: HagerZhang
            ls.mayterminate[] = false
        end
    else
        # Pick the initial step size according to HZ #I1-2
        if Tls <: HagerZhang
            mayterminate = ls.mayterminate
        else
            mayterminate = Ref{Bool}(false)
        end
        T = eltype(state.alpha)
        state.alpha = _hzI12(state.alpha, df, state.x, state.s, state.x_ls, phi_0, dphi_0,
                   is.ψ1, is.ψ2, is.ψ3, T(is.αmax), is.verbose, is.quadstep, mayterminate)
    end
    return state.alpha
end

# Pick the initial step size (HZ #I1-I2)
function _hzI12(alpha::T,
                df,
                x::AbstractArray{Tx},
                s::AbstractArray{Tx},
                x_new::AbstractArray{Tx},
                phi_0::T,
                dphi_0::T,
                psi1::Real,
                psi2::Real,
                psi3::Real,
                alphamax::T,
                verbose::Bool,
                quadstep::Bool,
                mayterminate) where {Tx,T}
    ϕ = make_ϕ(df, x_new, x, s)

    # Prevent values of `x_new` that are likely to make
    # ϕ(x_new) infinite
    iterfinitemax::Int = ceil(Int, -log2(eps(T)))

    alphatest = psi1 * alpha
    alphatest = min(alphatest, alphamax)
    alphatest == 0 && return alphatest
    phitest = ϕ(alphatest)
    alphatest, phitest = get_finite(alphatest, phitest, ϕ, psi3, iterfinitemax, phi_0, mayterminate)
    alphatest == 0 && return alphatest

    mayterminate[] = quadstep_success = false
    if quadstep
        a = ((phitest - phi_0) / alphatest - dphi_0) / alphatest  # quadratic fit
        if verbose == true
            println("quadfit: alphatest = ", alphatest,
                    ", phi_0 = ", phi_0,
                    ", dphi_0 = ", dphi_0,
                    ", phitest = ", phitest,
                    ", quadcoef = ", a)
        end
        if isfinite(a) && a > 0 && phitest <= phi_0
            alphatest2 = -dphi_0 / 2 / a # if convex, choose minimum of quadratic
            alphatest2 = min(alphatest2, alphamax)
            phitest2 = ϕ(alphatest2)
            if isfinite(phitest2)
                mayterminate[] = quadstep_success = true
                alphatest = alphatest2
                phitest = phitest2
                if verbose == true
                    println("alpha guess (quadratic): ", alphatest,
                            ",(mayterminate = ", mayterminate[], ")")
                end
            end
        end
    end
    if (!quadstep || !quadstep_success) && phitest <= phi_0
        # If no quadstep or it fails, expand the interval.
        # While the phitest <= phi_0 condition was not in the paper, it gives a significant boost to the speed. The rationale behind it is that since the slope at alpha = 0 is negative, if phitest > phi_0 then a local minimum must be between alpha = 0 and alpha = alphatest, so alpha_test is good enough to return.
        alphatest = psi2 * alpha 
        alphatest = min(alphatest, alphamax)
        phitest = ϕ(alphatest)
        alphatest, phitest = get_finite(alphatest, phitest, ϕ, psi3, iterfinitemax, phi_0, mayterminate)
        if verbose == true
            println("alpha guess (expand): ", alphatest)
        end
    end
    return alphatest
end

# Generate initial guess for step size (HZ, stage I0)
function _hzI0(x::AbstractArray{Tx},
               gr::AbstractArray{Tx},
               f_x::T,
               alphamax::T,
               psi0::T = convert(T, 1)/100) where {Tx,T}
    zeroT = convert(T, 0)
    alpha = convert(T, 1)
    gr_max = maximum(abs, gr)
    if gr_max != zeroT
        x_max = maximum(abs, x)
        if x_max != zeroT
            alpha = psi0 * x_max / gr_max
        elseif f_x != zeroT
            alpha = psi0 * abs(f_x) / norm(gr)^2
        end
    end
    return min(alpha, alphamax)
end

function get_finite(alpha::T, phi, ϕ, factor, itermax, phi_0, mayterminate) where {T}
    iter = 1
    while !isfinite(phi)
        if iter >= itermax
            mayterminate[] = true
            return T(0), phi_0
        end

        alpha = factor * alpha
        phi = ϕ(alpha)
        iter += 1
    end
    return alpha, phi
end
