"""
    BackTracking(; c_1=1e-4, ρ_hi=0.5, ρ_lo=0.1, iterations=1_000, order=3, maxstep=Inf, cache=nothing)

A backtracking line-search that uses a quadratic or cubic interpolant to
determine the reduction in step-size. E.g., if `f(α) > f(0) + c₁ α f'(0)`, then
the quadratic interpolant of `f(0)`, `f'(0)`, `f(α)` has a minimiser `α'` in the
open interval `(0, α)`. More strongly, there exists a factor `ρ = ρ(c_1)` such
that `α' ≦ ρ α`.

This is a modification of the algorithm described in Nocedal Wright (2nd ed),
Sec. 3.5.
"""
@kwdef struct BackTracking{TF, TI} <: AbstractLineSearch
    c_1::TF = 1e-4
    ρ_hi::TF = 0.5
    ρ_lo::TF = 0.1
    iterations::TI = 1_000
    order::TI = 3
    maxstep::TF = Inf
    cache::Union{Nothing,LineSearchCache{TF}} = nothing

    function BackTracking{TF,TI}(c_1, ρ_hi, ρ_lo, iterations, order, maxstep, cache) where {TF,TI}
        if !(0 < c_1 < 1)
            throw(ArgumentError(
                LazyString("The Armijo constant must satisfy 0 < c_1 < 1. Got c_1 = ", c_1),
            ))
        end
        if !(0 < ρ_lo <= ρ_hi < 1)
            throw(ArgumentError(
                LazyString("The backtracking factors must satisfy 0 < ρ_lo <= ρ_hi < 1. Got ρ_lo = ", ρ_lo, " and ρ_hi = ", ρ_hi),
            ))
        end
        if !(order in (2, 3))
            throw(ArgumentError(LazyString("The interpolation order must be 2 or 3. Got order = ", order)))
        end
        if !(iterations > 0)
            throw(ArgumentError(LazyString("The iteration limit must be positive. Got iterations = ", iterations)))
        end
        if !(maxstep > 0)
            throw(ArgumentError(LazyString("The maximum step length must be positive. Got maxstep = ", maxstep)))
        end
        return new{TF,TI}(c_1, ρ_hi, ρ_lo, iterations, order, maxstep, cache)
    end
end
BackTracking(c_1::TF, ρ_hi::TF, ρ_lo::TF, iterations::TI, order::TI, maxstep::TF,
             cache::Union{Nothing,LineSearchCache{TF}}) where {TF,TI} =
    BackTracking{TF,TI}(c_1, ρ_hi, ρ_lo, iterations, order, maxstep, cache)
BackTracking{TF}(args...; kwargs...) where TF = BackTracking{TF,Int}(args...; kwargs...)

function (ls::BackTracking)(df::AbstractObjective, x::AbstractArray{T}, s::AbstractArray{T},
                            α_0::Tα = real(T)(1), x_new::AbstractArray{T} = similar(x), ϕ_0 = nothing, dϕ_0 = nothing, alphamax = typemax(real(T))) where {T, Tα}
    # The search itself only ever needs ϕ; the others fill in a missing endpoint
    ϕ = make_ϕ(df, x_new, x, s)

    if isnothing(ϕ_0) && isnothing(dϕ_0)
        ϕ_0, dϕ_0 = make_ϕdϕ(df, x_new, x, s)(zero(Tα))
    elseif isnothing(ϕ_0)
        ϕ_0 = ϕ(zero(Tα))
    elseif isnothing(dϕ_0)
        dϕ_0 = make_dϕ(df, x_new, x, s)(zero(Tα))
    end

    α_0 = min(α_0, min(alphamax, ls.maxstep / norm(s, Inf)))
    ls(ϕ, α_0, ϕ_0, dϕ_0)
end

(ls::BackTracking)(ϕ, dϕ, ϕdϕ, αinitial, ϕ_0, dϕ_0) = ls(ϕ, αinitial, ϕ_0, dϕ_0)

# TODO: Should we deprecate the interface that only uses the ϕ argument?
function (ls::BackTracking)(ϕ, αinitial::Tα, ϕ_0, dϕ_0) where Tα
    (; c_1, ρ_hi, ρ_lo, iterations, order, cache) = ls
    emptycache!(cache)
    pushcache!(cache, 0, ϕ_0, dϕ_0)  # backtracking doesn't use the slope except here

    if !(isfinite(ϕ_0) && isfinite(dϕ_0))
        throw(LineSearchException("Value and slope at step length = 0 must be finite.", zero(Tα)))
    end
    # dϕ_0 == 0 is accepted: Armijo then reduces to a plain decrease test
    if dϕ_0 > 0
        throw(LineSearchException("Search direction is not a direction of descent.", zero(Tα)))
    end

    iterfinitemax = -log2(eps(real(Tα)))

    # Count the total number of iterations
    iteration = 0

    ϕx_0, ϕx_1 = ϕ_0, ϕ_0

    α_1, α_2 = αinitial, αinitial

    ϕx_1 = ϕ(α_1)

    # Hard-coded backtrack until we find a finite function value
    iterfinite = 0
    while !isfinite(ϕx_1) && iterfinite < iterfinitemax
        iterfinite += 1
        α_1 = α_2
        α_2 = α_1/2

        ϕx_1 = ϕ(α_2)
    end
    pushcache!(cache, α_2, ϕx_1)
    if !isfinite(ϕx_1)
        throw(LineSearchException("Backtracking: failed to achieve finite new evaluation point.", α_2))
    end

    # Backtrack until we satisfy sufficient decrease condition
    while ϕx_1 > ϕ_0 + c_1 * α_2 * dϕ_0
        # Increment the number of steps we've had to perform
        iteration += 1

        # Ensure termination
        if iteration > iterations
            throw(LineSearchException("Linesearch failed to converge, reached maximum iterations $(iterations).",
                                      α_2))
        end

        # Shrink proposed step-size:
        if order == 2 || iteration == 1
            # backtracking via quadratic interpolation:
            # This interpolates the available data
            #    f(0), f'(0), f(α)
            # with a quadractic which is then minimised; this comes with a
            # guaranteed backtracking factor 0.5 * (1-c_1)^{-1} which is < 1
            # provided that c_1 < 1/2; the backtrack_condition at the beginning
            # of the function guarantees at least a backtracking factor ρ.
            α_tmp = - (dϕ_0 * α_2^2) / ( 2 * (ϕx_1 - ϕ_0 - dϕ_0*α_2) )
        else
            div = one(Tα) / (α_1^2 * α_2^2 * (α_2 - α_1))
            a = (α_1^2*(ϕx_1 - ϕ_0 - dϕ_0*α_2) - α_2^2*(ϕx_0 - ϕ_0 - dϕ_0*α_1))*div
            b = (-α_1^3*(ϕx_1 - ϕ_0 - dϕ_0*α_2) + α_2^3*(ϕx_0 - ϕ_0 - dϕ_0*α_1))*div

            if isapprox(a, zero(a), atol=eps(real(Tα)))
                α_tmp = -dϕ_0 / (2*b)
            else
                # discriminant
                d = max(b^2 - 3*a*dϕ_0, zero(Tα))
                # quadratic equation root, rewritten via the conjugate to avoid
                # catastrophic cancellation when -b and sqrt(d) are close
                α_tmp = -dϕ_0 / (b + sqrt(d))
            end
        end

        α_1 = α_2

        α_tmp = NaNMath.min(α_tmp, α_2*ρ_hi) # avoid too small reductions
        α_2 = NaNMath.max(α_tmp, α_2*ρ_lo) # avoid too big reductions

        # Evaluate f(x) at proposed position. If non-finite, halve α
        # until we recover a finite value. This really shouldn't happen,
        # but some objective functions with ODE solves etc inside can somehow
        # unexpectedly fail next to fine candidates.
        ϕx_new = ϕ(α_2)
        pushcache!(cache, α_2, ϕx_new)
        while !isfinite(ϕx_new) && iteration < iterations
            iteration += 1
            α_2 = α_2 / 2
            ϕx_new = ϕ(α_2)
            pushcache!(cache, α_2, ϕx_new)
        end
        if !isfinite(ϕx_new)
            throw(LineSearchException("Backtracking: failed to achieve finite new evaluation point.", α_2))
        end
        ϕx_0, ϕx_1 = ϕx_1, ϕx_new
    end

    return α_2, ϕx_1
end
