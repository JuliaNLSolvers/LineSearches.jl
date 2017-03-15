"""
`backtracking!` is a backtracking line-search that uses
a quadratic or cubic interpolant to determine the reduction in step-size.
E.g.,
if f(α) > f(0) + c₁ α f'(0), then the quadratic interpolant of
f(0), f'(0), f(α) has a minimiser α' in the open interval (0, α). More strongly,
there exists a factor ρ = ρ(c₁) such that α' ≦ ρ α.

This is a modification of the algorithm described in Nocedal Wright (2nd ed), Sec. 3.5.
"""

bt3!{T}(df,
        x::Vector{T},
        s::Vector,
        x_scratch::Vector,
        gr_scratch::Vector,
        lsr::LineSearchResults,
        alpha::Real = 1.0,
        mayterminate::Bool = false,
        c1::Real = 1e-4,
        rhohi::Real = 0.5,
        rholo::Real = 0.1,
        iterations::Integer = 1_000) =
            backtracking!(df,x,s,x_scratch,gr_scratch,
                          lsr,alpha,mayterminate,c1,
                          rhohi,rholo,iterations,3)
bt2!{T}(df,
        x::Vector{T},
        s::Vector,
        x_scratch::Vector,
        gr_scratch::Vector,
        lsr::LineSearchResults,
        alpha::Real = 1.0,
        mayterminate::Bool = false,
        c1::Real = 1e-4,
        rhohi::Real = 0.5,
        rholo::Real = 0.1,
        iterations::Integer = 1_000) =
            backtracking!(df,x,s,x_scratch,gr_scratch,
                          lsr,alpha,mayterminate,c1,
                          rhohi,rholo,iterations,2)


function backtracking!{T}(df,
                          x::Vector{T},
                          s::Vector,
                          x_scratch::Vector,
                          gr_scratch::Vector,
                          lsr::LineSearchResults,
                          alpha::Real = 1.0,
                          mayterminate::Bool = false,
                          c1::Real = 1e-4,
                          rhohi::Real = 0.5,
                          rholo::Real = 0.1,
                          iterations::Integer = 1_000,
                          order::Int = 3)
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

    # Track calls to function and gradient
    f_calls = 0
    g_calls = 0

    # Count number of parameters
    n = length(x)

    # read f_x and slope from LineSearchResults
    f_x = lsr.value[end]
    gxp = lsr.slope[end]

    # Tentatively move a distance of alpha in the direction of s
    @simd for i in 1:n
        @inbounds x_scratch[i] = x[i] + alpha * s[i]
    end
    push!(lsr.alpha, alpha)

    # Backtrack until we satisfy sufficient decrease condition
    f_x_scratch = NLSolversBase.value!(df, x_scratch)
    push!(lsr.value, f_x_scratch)
    f_calls += 1
    while f_x_scratch > f_x + c1 * alpha * gxp
        # Increment the number of steps we've had to perform
        iteration += 1

        # Ensure termination
        if iteration > iterations
            throw(LineSearchException("Linesearch failed to converge, reached maximum iterations $(iterations).",
                                      lsr.alpha[end], f_calls, g_calls,lsr))
        end

        # Shrink proposed step-size:
        if order == 2 || iteration == 1
            # backtracking via interpolation:
            # This interpolates the available data
            #    f(0), f'(0), f(α)
            # with a quadractic which is then minimised; this comes with a
            # guaranteed backtracking factor 0.5 * (1-c1)^{-1} which is < 1
            # provided that c1 < 1/2; the backtrack_condition at the beginning
            # of the function guarantees at least a backtracking factor rho.
            alphatmp = - (gxp * alpha^2) / ( 2.0 * (f_x_scratch - f_x - gxp*alpha) )
        else
            alpha0 = lsr.alpha[end-1]
            alpha1 = lsr.alpha[end]
            phi0 = lsr.value[end-1]
            phi1 = lsr.value[end]

            div = 1.0/(alpha0^2*alpha1^2*(alpha1-alpha0))
            a = (alpha0^2*(phi1-f_x-gxp*alpha1)-alpha1^2*(phi0-f_x-gxp*alpha0))*div
            b = (-alpha0^3*(phi1-f_x-gxp*alpha1)+alpha1^3*(phi0-f_x-gxp*alpha0))*div

            if isapprox(a,0)
                alphatmp = gxp / (2.0*b)
            else
                discr = max(b^2-3*a*gxp, 0.)
                alphatmp = (-b + sqrt(discr)) / (3.0*a)
            end
        end
        alphatmp =  min(alphatmp, alpha*rhohi) # avoid too small reductions
        alpha = max(alphatmp, alpha*rholo) # avoid too big reductions

        push!(lsr.alpha, alpha)

        # Update proposed position
        @simd for i in 1:n
            @inbounds x_scratch[i] = x[i] + alpha * s[i]
        end

        # Evaluate f(x) at proposed position
        f_x_scratch = NLSolversBase.value!(df, x_scratch)
        f_calls += 1
        push!(lsr.value, f_x_scratch)
    end

    return alpha, f_calls, g_calls
end
