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

function (is::InitialStatic{T})(state, dphi0, df) where T
    state.alpha = is.alpha
    if is.scaled == true && (ns = norm(state.s)) > zero(T)
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

function (is::InitialPrevious)(state, dphi0, df)
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

The proposed step length lies between `αmin` and 1.0, where
`αmin = 0.1` per default. (This value is somewhat arbitrary).

Sets initial guess to the parameter `α0` at the first iteration (default 1.0).
"""
@with_kw struct InitialQuadratic{T}
    αmin::T = 0.1 # Minimum initial step size (value somewhat arbitrary)
    α0::T   = 1.0 # Fallback at first iteration
end

function (is::InitialQuadratic{T})(state, dphi0, df) where T
    if isnan(state.f_x_previous) || dphi0 ≈ zero(T)
        # If we're at the first iteration
        αguess = is.α0
    else
        αguess = 2.0 * (NLSolversBase.value(df) - state.f_x_previous) / dphi0
        αguess = min(one(T), 1.01*αguess) # See Nocedal+Wright
        αguess = max(is.αmin, αguess)
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

The proposed step length lies between `αmin` and 1.0, where
`αmin = 0.1` per default. (This value is somewhat arbitrary).

Sets initial guess to the parameter `α0` at the first iteration (default 1.0).
"""
@with_kw struct InitialConstantChange{T}
    αmin::T = 0.1 # Minimum initial step size (value somewhat arbitrary)
    α0::T   = 1.0 # Fallback at first iteration
end

function (is::InitialConstantChange{T})(state, dphi0, df) where T
    if isnan(state.dphi0_previous) || dphi0 ≈ zero(T)
        # If we're at the first iteration
        αguess = is.α0
    else
        # state.alpha is the previously used step length
        αguess = state.alpha * state.dphi0_previous / dphi0
        αguess = min(one(T), 1.01*αguess) # See Nocedal+Wright
        αguess = max(is.αmin, αguess)
    end
    state.alpha = αguess
    # TODO: Should `mayterminate` be true or false? Does it depend on which line search we use?
    state.mayterminate = false
end
