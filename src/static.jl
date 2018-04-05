# This linesearch does not perform a search at all, but takes the step
# alpha in the search direction. By default, alpha = 1.0 is the full step.
# This algorithm is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.
# NOTE: alpha passed from outside is ignored here, but can be passed to allow
# for generic use. We use ls.alpha instead.

"""
`Static`: defines a static linesearch, i.e. with fixed step-size. E.g., initialise
with `Static(alpha = 0.3141)` for fixed step-size 0.3141. Default is 1.0.
"""
@with_kw struct Static{T}
    alpha::T = 1.0
    scaled::Bool = false # Scales step. alpha ← min(alpha,||s||_2) / ||s||_2
end

function (ls::Static)(df::AbstractObjective, x, s, α, x_new = similar(x), ϕ_0 = nothing, dϕ_0 = nothing)
    ϕ = make_ϕ(df, x_new, x, s)
    ls(ϕ, α)
end

function (ls::Static)(ϕ, alpha = 1.0)
    @unpack alpha, scaled = ls
    @assert alpha > 0 # This should really be done at the constructor level

    ϕα = ϕ(alpha)

    return alpha, ϕα
end
