"""
    Static()

A static linesearch that accepts the initial step length unchanged (no parameters
to configure).

`Static` is intended for methods with well-scaled updates; i.e. Newton, on
well-behaved problems.
"""
struct Static <: AbstractLineSearch end

function (ls::Static)(df::AbstractObjective, x, s, α, x_new = similar(x), ϕ_0 = nothing, dϕ_0 = nothing)
    ϕ = make_ϕ(df, x_new, x, s)
    α_star, ϕα = ls(ϕ, α)
    set_x_new!(x_new, x, s, α_star)
    return α_star, ϕα
end

(ls::Static)(ϕ, dϕ, ϕdϕ, α, ϕ_0, dϕ_0) = ls(ϕ, α)

# TODO: Should we deprecate the interface that only uses the ϕ argument?
function (ls::Static)(ϕ, α::Tα) where Tα
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
