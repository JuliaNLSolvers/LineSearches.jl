# This linesearch does not perform a search at all, but takes the step
# alpha in the search direction. By default, alpha = 1.0 is the full step.
# This algorithm is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.
# NOTE: alpha passed from outside is ignored here, but can be passed to allow
# for generic use. We use ls.alpha instead.

"""
`Static`: defines a static linesearch, i.e. with fixed step-size. E.g., initialise
with `Static(alpha = 0.3141)` for fixed step-size 0.3141. Default is 1.0.

You can also make this independent of the size of the step `s`, by using
`Static(scaled = true)`.
This will then use a steps-size alpha ← min(alpha,||s||_2) / ||s||_2
"""
@with_kw struct Static{T}
    alpha::T = 1.0
    scaled::Bool = false # Scales step. alpha ← min(alpha,||s||_2) / ||s||_2
end

(ls::Static)(df, x, s, x_new, phi0, dphi0, alpha, mayterminate) = (ls::Static)(df, x, s, x_new)

function (ls::Static)(df, x, s, x_new)
    @unpack alpha, scaled = ls
    @assert alpha > 0 # This should really be done at the constructor level

    ϕ = make_ϕ(df, x_new, x, s)

    if scaled == true && (ns = vecnorm(s)) > zero(typeof(alpha))
        alpha = min(alpha, ns) / ns
    end

    # All line searches are assume to have evaluated the function at
    # the last value. We will change this to actually be returned instead
    ϕx = ϕ(alpha)

    return alpha
end
