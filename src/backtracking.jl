"""
`BackTracking` specifies a backtracking line-search that uses
a quadratic or cubic interpolant to determine the reduction in step-size.
E.g.,
if f(α) > f(0) + c₁ α f'(0), then the quadratic interpolant of
f(0), f'(0), f(α) has a minimiser α' in the open interval (0, α). More strongly,
there exists a factor ρ = ρ(c₁) such that α' ≦ ρ α.

This is a modification of the algorithm described in Nocedal Wright (2nd ed), Sec. 3.5.
"""
@with_kw struct BackTracking{TF, TI}
    c_1::TF = 1e-4
    ρ_hi::TF = 0.5
    ρ_lo::TF = 0.1
    iterations::TI = 1_000
    order::TI = 3
    maxstep::TF = Inf
end

function (ls::BackTracking)(df::AbstractObjective, x::AbstractArray{T}, s::AbstractArray{T},
                            α_0::Tα = T(1), x_new::AbstractArray{T} = similar(x), ϕ_0 = nothing, dϕ_0 = nothing, alphamax = convert(T, Inf)) where {T, Tα}
    ϕ, dϕ = make_ϕ_dϕ(df, x_new, x, s)

    if ϕ_0 == nothing
        ϕ_0 = ϕ(α_0)
    end
    if dϕ_0 == nothing
        dϕ_0 = ϕ(α_0)
    end

    ls(ϕ, x, s, α_0, ϕ_0, dϕ_0, alphamax)
end
function (ls::BackTracking)(ϕ, x::AbstractArray{T}, s::AbstractArray{T}, α_0::Tα,
                            ϕ_0, dϕ_0, alphamax = convert(T, Inf)) where {T, Tα}

    @unpack c_1, ρ_hi, ρ_lo, iterations, order, maxstep = ls

    iterfinitemax = -log2(eps(T))

    @assert order in (2,3)
    # Check the input is valid, and modify otherwise
    #backtrack_condition = 1.0 - 1.0/(2*ρ) # want guaranteed backtrack factor
    #if c_1 >= backtrack_condition
    #    warn("""The Armijo constant c_1 is too large; replacing it with
    #                   $(backtrack_condition)""")
    #   c_1 = backtrack_condition
    #end

    # Count the total number of iterations
    iteration = 0

    ϕx_0, ϕx_1 = ϕ_0, ϕ_0

    α_1, α_2 = α_0, α_0

    ϕx_1 = ϕ(α_1)

    iterfinite = 0

    # Backtrack until we satisfy sufficient decrease condition
    while !isfinite(ϕx_1) && iterfinite < iterfinitemax
        iterfinite += 1
        α_1 = α_2
        α_2 = (T(1)/2)*α_1

        # Backtrack until we satisfy sufficient decrease condition
        ϕx_1 = ϕ(α_2)
    end

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

            if isapprox(a, zero(a), atol=eps(T))
                α_tmp = dϕ_0 / (2*b)
            else
                # discriminant
                d = max(b^2 - 3*a*dϕ_0, zero(Tα))
                # quadratic equation root
                α_tmp = (-b + sqrt(d)) / (3*a)
            end
        end
        α_tmp = NaNMath.min(α_tmp, α_2*ρ_hi) # avoid too small reductions
        α_tmp = NaNMath.max(α_tmp, α_2*ρ_lo) # avoid too big reductions

        # enforce a maximum step alpha * s (application specific, default is Inf)
        α_1 = α_2
        α_2 = min(α_tmp, maxstep / vecnorm(s, Inf))

        # Evaluate f(x) at proposed position
        ϕx_0, ϕx_1 = ϕx_1, ϕ(α_2)
    end

    return α_2, ϕx_1
end
