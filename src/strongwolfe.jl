# TODO: Implement safeguards


"""
`StrongWolfe`: This linesearch algorithm guarantees that the step length
satisfies the (strong) Wolfe conditions.
See Nocedal and Wright - Algorithms 3.5 and 3.6

This algorithm is mostly of theoretical interest, users should most likely
use `MoreThuente`, `HagerZhang` or `BackTracking`.

## Parameters:  (and defaults)
* `c1 = 1e-4`: Armijo condition
* `c2 = 0.9` : second (strong) Wolfe condition
* `rho = 2.0` : bracket growth
"""
@with_kw immutable StrongWolfe{T}
   c1::T = 1e-4
   c2::T = 0.9
   rho::T = 2.0
end

(ls::StrongWolfe)(args...) =
        strongwolfe!(args...; c1=ls.c1, c2=ls.c2, rho=ls.rho)

function strongwolfe!{T}(df,
                         x::Vector,
                         p::Vector,
                         x_new::Vector,
                         gr_new::Vector,
                         lsr::LineSearchResults{T},
                         alpha0::Real,
                         mayterminate::Bool;
                         c1::Real = 1e-4,
                         c2::Real = 0.9,
                         rho::Real = 2.0)
    # TODO: do we need gr_new anymore?
    # Any call to gradient! or value_gradient! would update df.g anyway

    # Parameter space
    n = length(x)

    # Step-sizes
    a_0 = 0.0
    a_iminus1 = a_0
    a_i = alpha0
    a_max = 65536.0

    # phi(alpha) = df.f(x + alpha * p)
    phi_0 = lsr.value[end]
    phi_a_iminus1 = phi_0
    phi_a_i = NaN

    # phi'(alpha) = vecdot(g(x + alpha * p), p)
    phiprime_0 = lsr.slope[end]
    phiprime_a_i = NaN

    # Iteration counter
    i = 1

    while a_i < a_max
        # Update x_new
        @simd for index in 1:n
            @inbounds x_new[index] = x[index] + a_i * p[index]
        end

        # Evaluate phi(a_i)
        phi_a_i = NLSolversBase.value!(df, x_new)

        # Test Wolfe conditions
        if (phi_a_i > phi_0 + c1 * a_i * phiprime_0) ||
            (phi_a_i >= phi_a_iminus1 && i > 1)
            a_star = zoom(a_iminus1, a_i,
                          phiprime_0, phi_0,
                          df, x, p, x_new, gr_new)
            return a_star
        end

        # Evaluate phi'(a_i)
        NLSolversBase.gradient!(df, x_new)
        gr_new[:] = gradient(df)

        phiprime_a_i = vecdot(gr_new, p)

        # Check condition 2
        if abs(phiprime_a_i) <= -c2 * phiprime_0
            return a_i
        end

        # Check condition 3
        if phiprime_a_i >= 0.0
            a_star = zoom(a_i, a_iminus1,
                          phiprime_0, phi_0,
                          df, x, p, x_new, gr_new)
            return a_star
        end

        # Choose a_iplus1 from the interval (a_i, a_max)
        a_iminus1 = a_i
        a_i *= rho

        # Update phi_a_iminus1
        phi_a_iminus1 = phi_a_i

        # Update iteration count
        i += 1
    end

    # Quasi-error response
    return a_max
end

function zoom(a_lo::Real,
              a_hi::Real,
              phiprime_0::Real,
              phi_0::Real,
              df,
              x::Vector,
              p::Vector,
              x_new::Vector,
              gr_new::Vector;
              c1::Real = 1e-4,
              c2::Real = 0.9)

    # Parameter space
    n = length(x)

    # Step-size
    a_j = NaN

    # Count iterations
    iteration = 0
    max_iterations = 10

    # Shrink bracket
    while iteration < max_iterations
        iteration += 1

        # Cache phi_a_lo
        @simd for index in 1:n
            @inbounds x_new[index] = x[index] + a_lo * p[index]
        end
        phi_a_lo = NLSolversBase.value_gradient!(df, x_new)
        gr_new[:] = NLSolversBase.gradient(df)
        phiprime_a_lo = vecdot(gr_new, p)

        # Cache phi_a_hi
        @simd for index in 1:n
            @inbounds x_new[index] = x[index] + a_hi * p[index]
        end
        phi_a_hi = NLSolversBase.value_gradient!(df, x_new)
        gr_new[:] = NLSolversBase.gradient(df)
        phiprime_a_hi = vecdot(gr_new, p)

        # Interpolate a_j
        if a_lo < a_hi
            a_j = interpolate(a_lo, a_hi,
                              phi_a_lo, phi_a_hi,
                              phiprime_a_lo, phiprime_a_hi)
        else
            # TODO: Check if this is needed
            a_j = interpolate(a_hi, a_lo,
                              phi_a_hi, phi_a_lo,
                              phiprime_a_hi, phiprime_a_lo)
        end

        # Update x_new
        @simd for index in 1:n
            @inbounds x_new[index] = x[index] + a_j * p[index]
        end

        # Evaluate phi(a_j)
        phi_a_j = NLSolversBase.value!(df, x_new)

        # Check Armijo
        if (phi_a_j > phi_0 + c1 * a_j * phiprime_0) ||
            (phi_a_j > phi_a_lo)
            a_hi = a_j
        else
            # Evaluate phiprime(a_j)
            NLSolversBase.gradient!(df, x_new)
            gr_new[:] = gradient(df)
            phiprime_a_j = vecdot(gr_new, p)

            if abs(phiprime_a_j) <= -c2 * phiprime_0
                return a_j
            end

            if phiprime_a_j * (a_hi - a_lo) >= 0.0
                a_hi = a_lo
            end

            a_lo = a_j
        end
    end

    # Quasi-error response
    return a_j
end

# a_lo = a_{i - 1}
# a_hi = a_{i}
function interpolate(a_i1::Real, a_i::Real,
                     phi_a_i1::Real, phi_a_i::Real,
                     phiprime_a_i1::Real, phiprime_a_i::Real)
    d1 = phiprime_a_i1 + phiprime_a_i -
        3.0 * (phi_a_i1 - phi_a_i) / (a_i1 - a_i)
    d2 = sqrt(d1 * d1 - phiprime_a_i1 * phiprime_a_i)
    return a_i - (a_i - a_i1) *
        ((phiprime_a_i + d2 - d1) /
         (phiprime_a_i - phiprime_a_i1 + 2.0 * d2))
end
