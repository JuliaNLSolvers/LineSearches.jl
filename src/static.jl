# TODO: Remove deprecation by v6.0.0
# - Delete `DeprecatedStatic` and the `Static(;kwargs...)` function
# - Replace string `NewStatic` with Static

##########################
## DELETE FROM THIS LINE
##########################
function Static(;kwargs...)
    if isempty(kwargs)
        return NewStatic()
    else
        Base.warn_once("`Static` no longer has any parameters.
Use together with `InitialStatic` if you wish to retain the old functionality.")
        return DeprecatedStatic(;kwargs...)
    end
end

"""
`DeprecatedStatic`: defines a static linesearch, i.e. with fixed step-size. E.g., initialise
with `DeprecatedStatic(alpha = 0.3141)` for fixed step-size 0.3141. Default is 1.0.
You can also make this independent of the size of the step `s`, by using
`DeprecatedStatic(scaled = true)`.
This will then use a step-size alpha ← min(alpha,||s||_2) / ||s||_2
"""
@with_kw struct DeprecatedStatic{T}
    alpha::T = 1.0
    scaled::Bool = false # Scales step. alpha ← min(alpha,||s||_2) / ||s||_2
end

function (ls::DeprecatedStatic)(df::AbstractObjective, x, s, α, x_new = similar(x), phi0 = nothing, dphi0 = nothing)
    ϕ = make_ϕ(df, x_new, x, s)
    ls(ϕ, α)
end

function (ls::DeprecatedStatic)(ϕ, α::Tα) where Tα
    @unpack α, scaled = ls
    @assert α > Tα(0) # This should really be done at the constructor level

    if scaled == true && (ns = vecnorm(s)) > Tα(0)
        α = min(α, ns) / ns
    end

    ϕα = ϕ(α)

    # Hard-coded backtrack until we find a finite function value
    iterfinite = 0
    iterfinitemax = -log2(eps(real(Tα)))
    while !isfinite(ϕα) && iterfinite < iterfinitemax
        iterfinite += 1
        αold = α
        α    = αold/2

        ϕα = ϕ(α)
    end

    return α, ϕα
end
##########################
## DELETE UNTIL THIS LINE
##########################

"""
`NewStatic`: defines a static linesearch which returns the initial step length.

`NewStatic` is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.
"""
immutable NewStatic end

function (ls::NewStatic)(df::AbstractObjective, x, s, α, x_new = similar(x), ϕ_0 = nothing, dϕ_0 = nothing)
    ϕ = make_ϕ(df, x_new, x, s)
    ls(ϕ, x, s, α)
end

function (ls::NewStatic)(ϕ, x, s, α::Tα) where Tα
    @assert α > real(Tα(0))
    ϕα = ϕ(α)

    # Hard-coded backtrack until we find a finite function value
    iterfinite = 0
    iterfinitemax = -log2(eps(real(Tα)))
    while !isfinite(ϕα) && iterfinite < iterfinitemax
        iterfinite += 1
        αold = α
        α    = αold/2

        ϕα = ϕ(α)
    end

    return α, ϕα
end
