# Code used to use phi, a 1-parameter function induced by f and s
# Need to pass s as an explicit parameter
function alphatry{T}(alpha::T,
                     df,
                     x::Array,
                     s::Array,
                     xtmp::Array,
                     gtmp::Array,
                     lsr::LineSearchResults,
                     psi1::Real = convert(T,0.2),
                     psi2::Real = convert(T,2),
                     psi3::Real = convert(T,0.1),
                     iterfinitemax::Integer = ceil(Integer, -log2(eps(T))),
                     alphamax::Real = convert(T, Inf),
                     verbose::Bool = false)

    phi0 = lsr.value[1]
    dphi0 = lsr.slope[1]

    alphatest = psi1 * alpha
    alphatest = min(alphatest, alphamax)

    # Use xtmp here
    phitest = NLSolversBase.value!(df, x + alphatest * s)

    iterfinite = 1
    while !isfinite(phitest)
        alphatest = psi3 * alphatest
        # Use xtmp here
        phitest = NLSolversBase.value!(df, x + alphatest * s)
        lsr.nfailures += 1
        iterfinite += 1
        if iterfinite >= iterfinitemax
            return zero(T), true
#             error("Failed to achieve finite test value; alphatest = ", alphatest)
        end
    end
    a = ((phitest-phi0)/alphatest - dphi0)/alphatest  # quadratic fit
    if verbose == true
        println("quadfit: alphatest = ", alphatest,
                ", phi0 = ", phi0,
                ", phitest = ", phitest,
                ", quadcoef = ", a)
    end
    mayterminate = false
    if isfinite(a) && a > 0 && phitest <= phi0
        alpha = -dphi0 / 2 / a # if convex, choose minimum of quadratic
        if alpha == 0
            error("alpha is zero. dphi0 = ", dphi0, ", phi0 = ", phi0, ", phitest = ", phitest, ", alphatest = ", alphatest, ", a = ", a)
        end
        if alpha <= alphamax
            mayterminate = true
        else
            alpha = alphamax
            mayterminate = false
        end
        if verbose == true
            println("alpha guess (quadratic): ", alpha,
                    ",(mayterminate = ", mayterminate, ")")
        end
    else
        if phitest > phi0
            alpha = alphatest
        else
            alpha *= psi2 # if not convex, expand the interval
        end
    end
    alpha = min(alphamax, alpha)
    if verbose == true
        println("alpha guess (expand): ", alpha)
    end
    return alpha, mayterminate
end


# Generate initial guess for step size (HZ, stage I0)
function alphainit{T}(alpha::Real,
                      x::Array{T},
                      gr::Array,
                      f_x::Real,
                      psi0::T = convert(T,0.01))
    if isnan(alpha)
        alpha = one(T)
        gr_max = maximum(abs(gr))
        if gr_max != 0.0
            x_max = maximum(abs(x))
            if x_max != 0.0
                alpha = psi0 * x_max / gr_max
            elseif f_x != 0.0
                alpha = psi0 * abs(f_x) / vecnorm(gr)
            end
        end
    end
    return alpha
end
