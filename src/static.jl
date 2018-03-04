# This linesearch does not perform a search at all, but takes the step
# alpha in the search direction. By default, alpha = 1.0 is the full step.
# This algorithm is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.

"""
`Static`: defines a static linesearch, i.e. with fixed step-size. E.g., initialise
with `Static(alpha = 0.3141)` for fixed step-size 0.3141. Default is 1.0.

You can also make this independent of the size of the step `s`, by using
`Static(scaled = true)`.
This will then use a steps-size alpha ← min(ls.alpha,||s||_2) / ||s||_2
"""
@with_kw struct Static{T}
    alpha::T = 1.0
    scaled::Bool = false # Scales step. alpha ← min(alpha,||s||_2) / ||s||_2
end

(ls::Static)(df, x, s, x_scratch, phi0, dphi0, alpha, mayterminate) = (ls::Static)(df, x, s, x_scratch, alpha)
function (ls::Static)(df, x, s, x_scratch, alpha)
    # NOTE: alpha is ignored here, and we use ls.alpha instead

    if ls.scaled == true && (ns = vecnorm(s)) > zero(typeof(ls.alpha))
        scaledalpha = min(ls.alpha, ns) / ns
        retval = _static!(df, x, s, x_scratch, scaledalpha)
    else
        retval = _static!(df, x, s, x_scratch, ls.alpha)
    end
    return retval
end


function _static!(df,
                x::AbstractArray{T},
                s::AbstractArray{T},
                x_scratch::AbstractArray{T},
                alpha::Real = 1.0) where T
    @assert alpha > 0

    # Count number of parameters
    n = length(x)
    # Move a distance of alpha in the direction of s
    x_scratch .= x .+ alpha.*s

    # Evaluate f(x) at new position
    f_x_scratch = NLSolversBase.value!(df, x_scratch)

    return alpha
end
