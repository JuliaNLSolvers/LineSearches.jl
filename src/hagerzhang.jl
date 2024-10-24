#
# Conjugate gradient line search implementation from:
#   W. W. Hager and H. Zhang (2006) Algorithm 851: CG_DESCENT, a
#     conjugate gradient method with guaranteed descent. ACM
#     Transactions on Mathematical Software 32: 113–137.
#
# Code comments such as "HZ, stage X" or "HZ, eqs Y" are with
# reference to a particular point in this paper.
#
# Enhancements:
#
# - checks for Inf/NaN function values
#
# - support maximum value of alpha (equivalently, c). This
#   facilitates using these routines for constrained minimization
#   when you can calculate the distance along the path to the
#   disallowed region. (When you can't easily calculate that
#   distance, it can still be handled by returning Inf/NaN for
#   exterior points. It's just more efficient if you know the
#   maximum, because you don't have to test values that won't
#   work.) The maximum should be specified as the largest value for
#   which a finite value will be returned.  See, e.g., limits_box
#   below.  The default value for alphamax is Inf. See alphamaxfunc
#   for cgdescent and alphamax for HagerZhang.

const DEFAULTDELTA = 0.1 # Values taken from HZ paper (Nocedal & Wright recommends 0.01?)
const DEFAULTSIGMA = 0.9 # Values taken from HZ paper (Nocedal & Wright recommends 0.1 for GradientDescent)

# NOTE:
#   [1] The type `T` in the `HagerZhang{T}` need not be the same `T` as in
#       `hagerzhang!{T}`; in the latter, `T` comes from the input vector `x`.
#   [2] the only method parameter that is not included in the
#       type is `iterfinitemax` since this value needs to be
#       inferred from the input vector `x` and not from the type information
#       on the parameters


"""
Conjugate gradient line search implementation from:
  W. W. Hager and H. Zhang (2006) Algorithm 851: CG_DESCENT, a
    conjugate gradient method with guaranteed descent. ACM
    Transactions on Mathematical Software 32: 113–137.
"""
@with_kw struct HagerZhang{T, Tm} <: AbstractLineSearch
   delta::T = DEFAULTDELTA # c_1 Wolfe sufficient decrease condition
   sigma::T = DEFAULTSIGMA # c_2 Wolfe curvature condition (Recommend 0.1 for GradientDescent)
   alphamax::T = Inf
   rho::T = 5.0
   epsilon::T = 1e-6
   epsilonk::Union{Nothing,Base.RefValue{T}} = nothing
   gamma::T = 0.66
   theta::T = 0.5
   linesearchmax::Int = 50
   psi3::T = 0.1
   display::Int = 0         # non-functional
   mayterminate::Tm = Ref{Bool}(false)
   cache::Union{Nothing,LineSearchCache{T}} = nothing
end
HagerZhang{T}(args...; kwargs...) where T = HagerZhang{T, Base.RefValue{Bool}}(args...; kwargs...)

function (ls::HagerZhang)(df::AbstractObjective, x::AbstractArray{T},
                          s::AbstractArray{T}, α::Real,
                          x_new::AbstractArray{T}, phi_0::Real, dphi_0::Real) where T
    ϕ, ϕdϕ = make_ϕ_ϕdϕ(df, x_new, x, s)
    ls(ϕ, ϕdϕ, α::Real, phi_0, dphi_0)
end

(ls::HagerZhang)(ϕ, dϕ, ϕdϕ, c, phi_0, dphi_0) = ls(ϕ, ϕdϕ, c, phi_0, dphi_0)
(ls::HagerZhang)(ϕ,     ϕdϕ, c, phi_0, dphi_0) = ls(ϕdϕ, c, phi_0, dphi_0)

# TODO: Should we deprecate the interface that only uses the ϕ and ϕdϕ arguments?
function (ls::HagerZhang)(ϕdϕ,
                          c::Tα,
                          phi_0::Tϕ,
                          dphi_0::Real) where {Tα, Tϕ}
    @unpack alphamax, gamma, linesearchmax, cache = ls
    if cache === nothing
        cache = LineSearchCache{Tα}()  # this won't be seen by the user, but it's used internally
    else
        eltype(cache) === Tα || throw(LineSearchException("Element type of line search cache doesn't match location type.", zero(Tα)))
        emptycache!(cache)
    end
    ϵₖ = Tα(ls.epsilonk !== nothing ? ls.epsilonk[] : ls.epsilon * abs(phi_0))
    γ = Tα(ls.gamma)

    pushcache!(cache, zero(Tα), phi_0, dphi_0)
    if !(isfinite(phi_0) && isfinite(dphi_0))
        throw(LineSearchException("Value and slope at step length = 0 must be finite.", zero(Tα)))
    end
    c > zero(Tα) || throw(LineSearchException("Line search step `c` must be positive.", zero(Tα)))
    if dphi_0 * eps(c) >= abs(phi_0)
        throw(LineSearchException("Search direction is not a direction of descent.", zero(Tα)))
    elseif dphi_0 >= zero(dphi_0)
        return zero(Tα), phi_0
    end

    j = 0
    aⱼ, bⱼ, terminate = bracket(ϕdϕ, c, phi_0, dphi_0, ϵₖ, ls, cache)           # L0
    while !terminate && j < linesearchmax
        a, b, terminate = secant²(ϕdϕ, aⱼ, bⱼ, phi_0, dphi_0, ϵₖ, ls, cache)    # L1
        terminate && break
        if b - a > γ * (bⱼ - aⱼ)                                                # L2
            c = (a + b) / 2
            a, b, terminate = update(ϕdϕ, a, b, c, phi_0, dphi_0, ϵₖ, ls, cache)
        end
        aⱼ, bⱼ = a, b                                                           # L3
        j += 1
    end
    j >= linesearchmax && throw(LineSearchException("Line search failed to converge, reached maximum iterations $(linesearchmax).", aⱼ))

    ls.mayterminate[] = false
    return cache.alphas[end], cache.values[end]
end

function satisfies_wolfe(ls::HagerZhang,
                         α::T,
                         ϕ_α::Real,
                         dϕ_α::Real,
                         ϕ_0::Real,
                         dϕ_0::Real,
                         ϵₖ::Real) where T<:Number
    δ, σ = ls.delta, ls.sigma
    δ * dϕ_0 * α >= ϕ_α - ϕ_0 && dϕ_α >= σ * dϕ_0 && return true            # T1
    return (2 * δ - 1) * dϕ_0 >= dϕ_α >= σ * dϕ_0 && ϕ_α <= ϕ_0 + ϵₖ        # T2
end

function bracket(ϕdϕ::Function,
                 c::Tα,
                 ϕ_0::Tϕ,
                 dϕ_0::Real,
                 ϵₖ::Real,
                 ls::HagerZhang,
                 cache::LineSearchCache) where {Tα, Tϕ}
    δ, σ, ρ = ls.delta, ls.sigma, ls.rho
    δ, σ, ρ = Tα(δ), Tα(σ), Tα(ρ)
    j, cⱼ = 0, c
    while true
        ϕ_j, dϕ_j = ϕdϕ(cⱼ)
        pushcache!(cache, cⱼ, ϕ_j, dϕ_j)
        if iszero(dϕ_j) || j > 0 || ls.mayterminate[]
            satisfies_wolfe(ls, cⱼ, ϕ_j, dϕ_j, ϕ_0, dϕ_0, ϵₖ) && return zero(Tα), cⱼ, true
        end
        j += 1
        if dϕ_j >= zero(dϕ_j)                                               # B1
            b = cⱼ
            i = j
            while true
                i == 0 && @show j cache.values ϕ_0 ϵₖ
                cache.values[i] < ϕ_0 + ϵₖ && break
                i -= 1
            end
            return cache.alphas[i], b, false
        end
        if ϕ_j > ϕ_0 + ϵₖ                                                   # B2
            return update3(ϕdϕ, zero(Tα), cⱼ, ϕ_0, dϕ_0, ϵₖ, ls, cache)
        end
        cⱼ = ρ * cⱼ
    end
end

# Step U3
function update3(ϕdϕ, ā, b̄, ϕ_0, dϕ_0, ϵₖ, ls, cache)
    θ = ls.theta
    while true
        d = (1 - θ) * ā + θ * b̄
        @assert ā < d < b̄
        ϕ_d, dϕ_d = ϕdϕ(d)
        pushcache!(cache, d, ϕ_d, dϕ_d)
        satisfies_wolfe(ls, d, ϕ_d, dϕ_d, ϕ_0, dϕ_0, ϵₖ) && return ā, d, true
        dϕ_d >= zero(dϕ_d) && return ā, d, false
        if ϕ_d <= ϕ_0 + ϵₖ
            ā = d
        else
            b̄ = d
        end
    end
end

function update(ϕdϕ, a, b, c, phi_0, dphi_0, ϵₖ, ls, cache)
    a < c < b || return a, b, false                                             # U0
    # Custom modification: sometimes c duplicates a previous value. Since
    # computing objectives and gradients can be expensive, it seems much better
    # to try to look it up.
    i = findfirst(==(c), cache.alphas)
    if i !== nothing
        ϕ_c, dϕ_c = cache.values[i], cache.slopes[i]
    else
        ϕ_c, dϕ_c = ϕdϕ(c)
        pushcache!(cache, c, ϕ_c, dϕ_c)
        satisfies_wolfe(ls, c, ϕ_c, dϕ_c, phi_0, dphi_0, ϵₖ) && return a, c, true
    end
    dϕ_c >= zero(dϕ_c) && return a, c, false                                    # U1
    ϕ_c <= phi_0 + ϵₖ && return c, b, false                                     # U2
    return update3(ϕdϕ, a, c, phi_0, dphi_0, ϵₖ, ls, cache)                     # U3
end

function secant²(ϕdϕ, a::T, b::T, phi_0, dphi_0, ϵₖ, ls, cache::LineSearchCache{T}) where T
    @assert a < b
    c = secant(a, b, cache)                                                    # S1
    A, B, terminate = update(ϕdϕ, a, b, c, phi_0, dphi_0, ϵₖ, ls, cache)
    terminate && return A, B, true
    c != A && c != B && return A, B, false                                     # S4, part b
    c̄ = c == B ? secant(b, B, cache) : secant(a, A, cache)                     # S2 & S3
    return update(ϕdϕ, A, B, c̄, phi_0, dphi_0, ϵₖ, ls, cache)                  # S4, part a
end

secant(a::T, b::T, cache::LineSearchCache{T}) where T =
    secant(a, b,
           cache.slopes[findfirst(==(a), cache.alphas)],
           cache.slopes[findfirst(==(b), cache.alphas)])

function secant(a::Real, b::Real, dϕ_a::Real, dϕ_b::Real)
    dϕ_a == dϕ_b && return (a + b)/2              # custom modification
    return (a * dϕ_b - b * dϕ_a) / (dϕ_b - dϕ_a)
end
