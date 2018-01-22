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
    c1::TF = 1e-4
    rhohi::TF = 0.5
    rholo::TF = 0.1
    iterations::TI = 1_000
    order::TI = 3
    maxstep::TF = Inf
end

(ls::BackTracking)(df, x, s, x_scratch, lsr, alpha, mayterminate) =
    _backtracking!(df, x, s, x_scratch, lsr, alpha, mayterminate,
             ls.c1, ls.rhohi, ls.rholo, ls.iterations, ls.order, ls.maxstep)


function _backtracking!(df,
                        x::AbstractArray{T},
                        s::AbstractArray{T},
                        x_scratch::AbstractArray{T},
                        lsr::LineSearchResults,
                        alpha::Real = 1.0,
                        mayterminate::Bool = false,
                        c1::Real = 1e-4,
                        rhohi::Real = 0.5,
                        rholo::Real = 0.1,
                        iterations::Integer = 1_000,
                        order::Int = 3,
                        maxstep::Real = Inf) where T
    iterfinitemax = -log2(eps(T))

    @assert order in (2,3)
    # Check the input is valid, and modify otherwise
    #backtrack_condition = 1.0 - 1.0/(2*rho) # want guaranteed backtrack factor
    #if c1 >= backtrack_condition
    #    warn("""The Armijo constant c1 is too large; replacing it with
    #                   $(backtrack_condition)""")
    #   c1 = backtrack_condition
    #end

    # Count the total number of iterations
    iteration = 0

    # Count number of parameters
    n = length(x)

    # read f_x and slope from LineSearchResults
    @inbounds f_x = lsr.value[end]
    @inbounds gxp = lsr.slope[end]

    # Tentatively move a distance of alpha in the direction of s
    x_scratch .= x .+ alpha.*s
    push!(lsr.alpha, alpha)

    # Backtrack until we satisfy sufficient decrease condition
    f_x_scratch = NLSolversBase.value!(df, x_scratch)
    push!(lsr.value, f_x_scratch)

    iterfinite = 0
    while !isfinite(f_x_scratch) && iterfinite < iterfinitemax
        iterfinite += 1
        alpha *= 0.5
        # Tentatively move a distance of alpha in the direction of s
        x_scratch .= x .+ alpha.*s
        push!(lsr.alpha, alpha)

        # Backtrack until we satisfy sufficient decrease condition
        f_x_scratch = NLSolversBase.value!(df, x_scratch)
        push!(lsr.value, f_x_scratch)
    end

    while f_x_scratch > f_x + c1 * alpha * gxp
        # Increment the number of steps we've had to perform
        iteration += 1

        # Ensure termination
        if iteration > iterations
            throw(LineSearchException("Linesearch failed to converge, reached maximum iterations $(iterations).",
                                      lsr.alpha[end], lsr))
        end

        # Shrink proposed step-size:
        if order == 2 || iteration == 1
            # backtracking via quadratic interpolation:
            # This interpolates the available data
            #    f(0), f'(0), f(α)
            # with a quadractic which is then minimised; this comes with a
            # guaranteed backtracking factor 0.5 * (1-c1)^{-1} which is < 1
            # provided that c1 < 1/2; the backtrack_condition at the beginning
            # of the function guarantees at least a backtracking factor rho.
            alphatmp = - (gxp * alpha^2) / ( 2.0 * (f_x_scratch - f_x - gxp*alpha) )
        else
            # Backtracking via cubic interpolation
            @inbounds alpha0 = lsr.alpha[end-1]
            @inbounds alpha1 = lsr.alpha[end]
            @inbounds phi0 = lsr.value[end-1]
            @inbounds phi1 = lsr.value[end]

            div = one(alpha) / (alpha0^2 * alpha1^2 * (alpha1 - alpha0))
            a = (alpha0^2*(phi1-f_x-gxp*alpha1)-alpha1^2*(phi0-f_x-gxp*alpha0))*div
            b = (-alpha0^3*(phi1-f_x-gxp*alpha1)+alpha1^3*(phi0-f_x-gxp*alpha0))*div

            if isapprox(a, zero(a))
                alphatmp = gxp / (2.0*b)
            else
                discr = max(b^2-3*a*gxp, zero(alpha))
                alphatmp = (-b + sqrt(discr)) / (3.0*a)
            end
        end
        alphatmp =  NaNMath.min(alphatmp, alpha*rhohi) # avoid too small reductions
        alphatmp = NaNMath.max(alphatmp, alpha*rholo) # avoid too big reductions

        # enforce a maximum step alpha * s (application specific, default is Inf)
        alpha = min(alphatmp, maxstep / vecnorm(s, Inf))

        push!(lsr.alpha, alpha)

        # Update proposed position
        x_scratch .= x .+ alpha.*s

        # Evaluate f(x) at proposed position
        f_x_scratch = NLSolversBase.value!(df, x_scratch)
        push!(lsr.value, f_x_scratch)
    end

    return alpha
end
