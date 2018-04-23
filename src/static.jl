# This linesearch does not perform a search at all, but takes the step
# alpha in the search direction.
# This algorithm is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.

"""
`Static`: defines a static linesearch which returns the initial step length.
"""
immutable Static end

function (ls::Static)(df::AbstractObjective, x, s, α, x_new = similar(x), phi0 = nothing, dphi0 = nothing)
    ϕ = make_ϕ(df, x_new, x, s)
    ls(ϕ, x, s, α)
end

function (ls::Static)(ϕ, x, s, alpha)
    ϕα = ϕ(alpha)

    return alpha, ϕα
end
