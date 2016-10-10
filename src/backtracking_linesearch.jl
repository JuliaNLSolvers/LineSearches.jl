
"""
`interpbacktrack_linesearch!` is a backtracking line-search that uses
a quadratic interpolant to determine the reduction in step-size. Specifically,
if f(α) > f(0) + c₁ α f'(0), then the quadratic interpolant of
f(0), f'(0), f(α) has a minimiser α' in the open interval (0, α). More strongly,
there exists a factor ρ = ρ(c₁) such that α' ≦ ρ α. This makes
the `interpbacktrack_linesearch!` a backtracking type linesearch.

This is a modification of the algorithm described in Nocedal Wright (2nd ed), Sec. 3.5.

**Default Parameters**

* `c1 = 0.2` : Armijo condition
* `rho = 0.9` : decrease `rho` to automatically adjust `c1` to obtain a stronger backtracking guarantee
* `mindecfact = 0.25` : specifies another safe-guard, α' ← max(α', mindecfact * α) to
   make sure not too small steps are taken.
"""
interpbacktrack_linesearch!{T}(df,
                               x::Vector{T},
                               s::Vector,
                               x_scratch::Vector,
                               gr_scratch::Vector,
                               lsr::LineSearchResults,
                               alpha::Real = 1.0,
                               mayterminate::Bool = false,
                               c1::Real = 0.2,
                               c2::Real = 0.9,
                               rho=0.9,
                               iterations::Integer = 1_000,
                               mindecfact=0.25) =
   backtracking_linesearch!(df,
                             x,
                             s,
                             x_scratch,
                             gr_scratch,
                             lsr,
                             alpha,
                             mayterminate,
                             c1,
                             c2,
                             rho,
                             iterations,
                             true,
                             mindecfact)


function backtracking_linesearch!{T}(df,
                                     x::Vector{T},
                                     s::Vector,
                                     x_scratch::Vector,
                                     gr_scratch::Vector,
                                     lsr::LineSearchResults,
                                     alpha::Real = 1.0,
                                     mayterminate::Bool = false,
                                     c1::Real = 1e-4,
                                     c2::Real = 0.9,
                                     rho::Real = 0.9,
                                     iterations::Integer = 1_000,
                                     interp::Bool = false,
                                     mindecfact=0.25)

    # Check the input is valid, and modify otherwise
    if interp   # this means we are coming from interpbacktrack_linesearch!
       backtrack_condition = 1.0 - 1.0/(2*rho) # want guaranteed backtrack factor
       if c1 >= backtrack_condition
           warn("""The Armijo constant c1 is too large; replacing it with
                   $(backtrack_condition)""")
           c1 = backtrack_condition
       end
       if rho <= mindecfact
           warn("""rho ($rho) <= mindecfact ($mindecfact); revert to standard
                  backtracking with constant factor rho = $rho""")
           interp = false
       end
    end

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

    # Backtrack until we satisfy sufficient decrease condition
    f_x_scratch = df.f(x_scratch)
    f_calls += 1
    while f_x_scratch > f_x + c1 * alpha * gxp
        # Increment the number of steps we've had to perform
        iteration += 1

        # Ensure termination
        if iteration > iterations
            error("Too many iterations in backtracking_linesearch!")
        end

        # Shrink proposed step-size:
        if !interp
           # standard backtracking:
           alpha *= rho
        else
           # backtracking via interpolation:
           # This interpolates the available data
           #    f(0), f'(0), f(α)
           # with a quadractic which is then minimised; this comes with a
           # guaranteed backtracking factor 0.5 * (1-c1)^{-1} which is < 1
           # provided that c1 < 1/2; the backtrack_condition at the beginning
           # of the function guarantees at least a backtracking factor rho.
           alpha1 = - (gxp * alpha) / ( 2.0 * ((f_x_scratch - f_x)/alpha - gxp) )
           alpha = max(alpha1, alpha * min(mindecfact, rho))  # avoid miniscule steps
        end

        # Update proposed position
        @simd for i in 1:n
            @inbounds x_scratch[i] = x[i] + alpha * s[i]
        end

        # Evaluate f(x) at proposed position
        f_x_scratch = df.f(x_scratch)
        f_calls += 1
    end

    return alpha, f_calls, g_calls
end
