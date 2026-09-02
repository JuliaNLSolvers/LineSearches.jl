# TODO: Implement safeguards


"""
    StrongWolfe(; c_1=1e-4, c_2=0.9, ρ=2.0, cache=nothing)

A line search algorithm that guarantees the step length satisfies the (strong)
Wolfe conditions. See Nocedal and Wright — Algorithms 3.5 and 3.6.

This algorithm is mostly of theoretical interest; users should most likely use
`MoreThuente`, `HagerZhang`, or `BackTracking`.

## Parameters (and defaults)

* `c_1 = 1e-4`: Armijo condition
* `c_2 = 0.9` : second (strong) Wolfe condition
* `ρ = 2.0`   : bracket growth
"""
@kwdef struct StrongWolfe{T} <: AbstractLineSearch
    c_1::T = 1e-4
    c_2::T = 0.9
    ρ::T = 2.0
    cache::Union{Nothing,LineSearchCache{T}} = nothing

    function StrongWolfe{T}(c_1, c_2, ρ, cache) where T
        if !(0 < c_1 < c_2 < 1)
            throw(ArgumentError(
                LazyString("The Wolfe constants must satisfy 0 < c_1 < c_2 < 1. Got c_1 = ", c_1, " and c_2 = ", c_2),
            ))
        end
        if !(ρ > 1)
            throw(ArgumentError(
                LazyString("The bracket growth factor must be larger than one, otherwise bracketing never terminates. Got ρ = ", ρ),
            ))
        end
        return new{T}(c_1, c_2, ρ, cache)
    end
end
StrongWolfe(c_1::T, c_2::T, ρ::T, cache::Union{Nothing,LineSearchCache{T}}) where T =
    StrongWolfe{T}(c_1, c_2, ρ, cache)

"""
    (ls::StrongWolfe)(df, x::AbstractArray, p::AbstractArray, alpha0::Real, x_new, ϕ_0, dϕ_0) -> alpha, ϕalpha

Given a differentiable function `df` (in the sense of `NLSolversBase.OnceDifferentiable` or
`NLSolversBase.TwiceDifferentiable`), a multidimensional starting point `x` and step `p`,
and a guess `alpha0` for the step length, find an `alpha` satisfying the strong Wolfe conditions.

See the one-dimensional method for additional details.
"""
function (ls::StrongWolfe)(df, x::AbstractArray{T},
                           p::AbstractArray{T}, α::Real, x_new::AbstractArray{T},
                           ϕ_0, dϕ_0) where T
    ϕdϕ = make_ϕdϕ(df, x_new, x, p)
    ls(ϕdϕ, α, ϕ_0, dϕ_0)
end

(ls::StrongWolfe)(ϕ, dϕ, ϕdϕ, alpha0, ϕ_0, dϕ_0) = ls(ϕdϕ, alpha0, ϕ_0, dϕ_0)

"""
    (ls::StrongWolfe)(ϕdϕ, alpha0, ϕ_0, dϕ_0) -> alpha, ϕalpha

Given a combined-evaluation function `ϕdϕ(alpha) -> (ϕ(alpha), dϕ(alpha))` and an
initial guess `alpha0`, identify a value of `alpha > 0` satisfying the strong Wolfe
conditions.

`ϕ_0` and `dϕ_0` are the value and derivative, respectively, of `ϕ` at `alpha = 0.`

Both `alpha` and `ϕ(alpha)` are returned.
"""
function (ls::StrongWolfe)(ϕdϕ, alpha0::T, ϕ_0, dϕ_0) where T<:Real
    (; c_1, c_2, ρ, cache) = ls
    emptycache!(cache)

    pushcache!(cache, zero(T), ϕ_0, dϕ_0)

    if !(isfinite(ϕ_0) && isfinite(dϕ_0))
        throw(LineSearchException("Value and slope at step length = 0 must be finite.", zero(T)))
    end
    # dϕ_0 == 0 is accepted: Armijo then reduces to a plain decrease test
    if dϕ_0 > zero(T)
        throw(LineSearchException("Search direction is not a direction of descent.", zero(T)))
    end

    # Step-sizes
    a_iminus1 = zero(T)
    a_i = alpha0
    a_max = convert(T, 65536)

    # ϕ(alpha) = df.f(x + alpha * p) and ϕ'(alpha) = dot(g(x + alpha * p), p) at
    # a_{i - 1}, which is the endpoint zoom brackets against
    ϕ_a_iminus1, dϕ_a_iminus1 = ϕ_0, dϕ_0

    # Iteration counter
    i = 1

    while a_i < a_max
        # Both zoom and the curvature test below need the slope, so fuse the two
        ϕ_a_i, dϕ_a_i = ϕdϕ(a_i)
        pushcache!(cache, a_i, ϕ_a_i, dϕ_a_i)

        # Test Wolfe conditions
        if (ϕ_a_i > ϕ_0 + c_1 * a_i * dϕ_0) ||
            (ϕ_a_i >= ϕ_a_iminus1 && i > 1)
            return zoom(a_iminus1, ϕ_a_iminus1, dϕ_a_iminus1,
                        a_i, ϕ_a_i, dϕ_a_i,
                        dϕ_0, ϕ_0, ϕdϕ, cache, c_1, c_2)
        end

        # Check condition 2
        if abs(dϕ_a_i) <= -c_2 * dϕ_0
            return a_i, ϕ_a_i
        end

        # Check condition 3
        if dϕ_a_i >= zero(T)
            return zoom(a_i, ϕ_a_i, dϕ_a_i,
                        a_iminus1, ϕ_a_iminus1, dϕ_a_iminus1,
                        dϕ_0, ϕ_0, ϕdϕ, cache, c_1, c_2)
        end

        # Choose a_iplus1 from the interval (a_i, a_max)
        a_iminus1, ϕ_a_iminus1, dϕ_a_iminus1 = a_i, ϕ_a_i, dϕ_a_i
        a_i *= ρ

        # Update iteration count
        i += 1
    end

    throw(LineSearchException("StrongWolfe: bracketing reached a_max=$a_max without satisfying Wolfe conditions.", a_max))
end

# The caller has evaluated both endpoints, and they are only ever reassigned to
# already-evaluated points, so no endpoint is ever recomputed here.
function zoom(a_lo::T, ϕ_a_lo::Real, ϕprime_a_lo::Real,
              a_hi::T, ϕ_a_hi::Real, ϕprime_a_hi::Real,
              dϕ_0::Real,
              ϕ_0::Real,
              ϕdϕ,
              cache,
              c_1::Real,
              c_2::Real) where T

    # Step-size
    a_j = convert(T, NaN)

    # Count iterations
    iteration = 0
    max_iterations = 10

    # Shrink bracket
    while iteration < max_iterations
        iteration += 1

        # Interpolate a_j
        if a_lo < a_hi
            a_j = interpolate(a_lo, a_hi,
                              ϕ_a_lo, ϕ_a_hi,
                              ϕprime_a_lo, ϕprime_a_hi)
        else
            # TODO: Check if this is needed
            a_j = interpolate(a_hi, a_lo,
                              ϕ_a_hi, ϕ_a_lo,
                              ϕprime_a_hi, ϕprime_a_lo)
        end

        # a_j becomes an endpoint in every branch below, so its slope is needed either way
        ϕ_a_j, ϕprime_a_j = ϕdϕ(a_j)
        pushcache!(cache, a_j, ϕ_a_j, ϕprime_a_j)

        # Check Armijo
        if (ϕ_a_j > ϕ_0 + c_1 * a_j * dϕ_0) ||
            (ϕ_a_j > ϕ_a_lo)
            a_hi, ϕ_a_hi, ϕprime_a_hi = a_j, ϕ_a_j, ϕprime_a_j
        else
            if abs(ϕprime_a_j) <= -c_2 * dϕ_0
                return a_j, ϕ_a_j
            end

            # Reads the bracket width before a_lo moves below
            if ϕprime_a_j * (a_hi - a_lo) >= zero(T)
                a_hi, ϕ_a_hi, ϕprime_a_hi = a_lo, ϕ_a_lo, ϕprime_a_lo
            end

            a_lo, ϕ_a_lo, ϕprime_a_lo = a_j, ϕ_a_j, ϕprime_a_j
        end
    end

    throw(LineSearchException("StrongWolfe: zoom reached maximum iterations ($max_iterations) without satisfying Wolfe conditions.", a_j))
end

# a_lo = a_{i - 1}
# a_hi = a_{i}
function interpolate(a_i1::Real, a_i::Real,
                     ϕ_a_i1::Real, ϕ_a_i::Real,
                     dϕ_a_i1::Real, dϕ_a_i::Real)
    d1, d2 = cubic_d1_d2(a_i1, ϕ_a_i1, dϕ_a_i1, a_i, ϕ_a_i, dϕ_a_i)
    # Unlike the cubic in backtracking.jl, this numerator does not cancel, so
    # rewriting it through its conjugate would cost accuracy rather than gain it
    a_j = a_i - (a_i - a_i1) *
        ((dϕ_a_i + d2 - d1) /
         (dϕ_a_i - dϕ_a_i1 + 2 * d2))
    # A vanishing denominator or a degenerate cubic can leave the bracket
    lo, hi = minmax(a_i1, a_i)
    return lo < a_j < hi ? a_j : (a_i1 + a_i) / 2
end
