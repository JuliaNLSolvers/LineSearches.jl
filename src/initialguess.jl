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

function (is::InitialStatic{T})(state, phi_0, dphi_0, df) where T
    state.alpha = is.alpha
    if is.scaled == true && (ns = vecnorm(state.s)) > zero(T)
        # TODO: Type instability if there's a type mismatch between is.alpha and ns
        state.alpha *= min(is.alpha, ns) / ns
    end
    # TODO: Should `mayterminate` be true or false? Does it depend on which line search we use?
    state.mayterminate = false
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

function (is::InitialPrevious)(state, phi_0, dphi_0, df)
    if isnan(state.alpha)
        state.alpha = is.alpha
        # TODO: Should `mayterminate` be true or false? Does it depend on which line search we use?
        state.mayterminate = false
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

function (is::InitialQuadratic{T})(state, phi_0, dphi_0, df) where T
    if !isfinite(state.f_x_previous) || dphi_0 ≈ zero(T)
        # If we're at the first iteration
        αguess = is.α0
    else
        αguess = 2.0 * (NLSolversBase.value(df) - state.f_x_previous) / dphi_0
        αguess = NaNMath.max(is.αmin, state.alpha*is.ρ, αguess)
        αguess = NaNMath.min(is.αmax, αguess)
        # if αguess ≈ 1, then make it 1 (Newton-type behaviour)
        if is.snap2one[1] ≤ αguess ≤ is.snap2one[2]
            αguess = one(state.alpha)
        end
    end
    state.alpha = αguess
    # TODO: Should `mayterminate` be true or false? Does it depend on which line search we use?
    state.mayterminate = false
end

"""
Constant first-order change approximation to determine initial step length.

** This requires that the optimization algorithm stores dphi0 from the previous iteration **
(dphi0_previous = vecdot(∇f_{k-1}, s_{k-1}), where s is the step direction.

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
@with_kw struct InitialConstantChange{T}
    αmin::T = 1e-12 # Minimum initial step size (value somewhat arbitrary)
    αmax::T = 1.0   # Maximum initial step size
    α0::T   = 1.0   # Fallback at first iteration
    ρ::T    = 0.25  # maximum decrease from previous step
    snap2one::Tuple{T,T} = (0.75, Inf) # Set everything in this (closed) interval to 1.0
end

function (is::InitialConstantChange{T})(state, phi_0, dphi_0, df) where T
    if !isfinite(state.dphi_0_previous) || !isfinite(state.alpha) || dphi_0 ≈ zero(T)
        # If we're at the first iteration
        αguess = is.α0
    else
        # state.alpha is the previously used step length
        αguess = state.alpha * state.dphi_0_previous / dphi_0
        αguess = NaNMath.max(is.αmin, state.alpha*is.ρ, αguess)
        αguess = NaNMath.min(is.αmax, αguess)
        # if αguess ≈ 1, then make it 1 (Newton-type behaviour)
        if is.snap2one[1] ≤ αguess ≤ is.snap2one[2]
            αguess = one(state.alpha)
        end
    end
    state.alpha = αguess
    # TODO: Should `mayterminate` be true or false? Does it depend on which line search we use?
    state.mayterminate = false
end
