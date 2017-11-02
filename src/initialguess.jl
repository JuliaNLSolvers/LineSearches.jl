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

function (is::InitialStatic)(state, dphi0, df)
    state.alpha = is.alpha
    if is.scaled == true && (ns = norm(state.s)) > zero(typeof(is.alpha))
        # TODO: Type instability if there's a type mismatch between is.alpha and ns
        state.alpha *= min(is.alpha, ns) / ns
    end
    # TODO: Should `mayterminate` be true or false? Does it depend on which line search we use?
    state.mayterminate = true
end

"""
Use previous step length as initial guess.

If state.alpha is NaN, then return fallback value is.alpha
"""
@with_kw struct InitialPrevious{T}
    alpha::T = 1.0
end

function (is::InitialPrevious)(state, dphi0, df)
    if isnan(state.alpha)
        state.alpha = is.alpha
        # TODO: Should `mayterminate` be true or false? Does it depend on which line search we use?
        state.mayterminate = true
    end
end
