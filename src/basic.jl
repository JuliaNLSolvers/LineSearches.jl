# This linesearch does not perform a search at all, but takes the step
# alpha in the search direction. By default, alpha = 1.0 is the full step.
# This algorithm is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.

"""
`Static`: defines a static linesearch, i.e. with fixed step-size. E.g., initialise
with `Static(alpha = 0.3141)` for fixed step-size 0.3141. Default is 1.0.
"""
@with_kw immutable Static{T}
    alpha::T = 1.0
end

(ls::Static)(df, x, s, x_scratch, gr_scratch, lsr, alpha, mayterminate) =
        _static!(df, x, s, x_scratch, gr_scratch, lsr, ls.alpha, mayterminate)

function _static!{T}(df,
                   x::Vector{T},
                   s::Vector,
                   x_scratch::Vector,
                   gr_scratch::Vector,
                   lsr::LineSearchResults,
                   alpha::Real = 1.0,
                   mayterminate::Bool = false)
    @assert alpha > 0
    push!(lsr.alpha, alpha)

    # Count number of parameters
    n = length(x)
    # Move a distance of alpha in the direction of s
    @simd for i in 1:n
        @inbounds x_scratch[i] = x[i] + alpha * s[i]
    end

    # Evaluate f(x) at new position
    f_x_scratch = NLSolversBase.value!(df, x_scratch)
    push!(lsr.value, f_x_scratch)

    return alpha
end
