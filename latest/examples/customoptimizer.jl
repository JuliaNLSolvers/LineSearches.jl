# # Using LineSearches without Optim/NLsolve
#
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`customoptimizer.ipynb`](@__NBVIEWER_ROOT_URL__examples/generated/customoptimizer.ipynb)
#-
#
# This tutorial shows you how to use the line search algorithms in `LineSearches`
# for your own optimization algorithm that is not part of `Optim` or `NLsolve`.
#
# Say we have written a gradient descent optimization algorithm but would like to
# experiment with different line search algorithms.
# The algorithm is implemented as follows.

function gdoptimize(f, g!, fg!, x0::AbstractArray{T}, linesearch,
                    maxiter::Int = 10000,
                    g_rtol::T = sqrt(eps(T)), g_atol::T = eps(T)) where T <: Number
    x = copy(x0)
    gvec = similar(x)
    g!(gvec, x)
    fx = f(x)

    gnorm = norm(gvec)
    gtol = max(g_rtol*gnorm, g_atol)

    ## Univariate line search functions
    ϕ(α) = f(x .+ α.*s)
    function dϕ(α)
        g!(gvec, x .+ α.*s)
        return vecdot(gvec, s)
    end
    function ϕdϕ(α)
        phi = fg!(gvec, x .+ α.*s)
        dphi = vecdot(gvec, s)
        return (phi, dphi)
    end

    s = similar(gvec) # Step direction

    iter = 0
    while iter < maxiter && gnorm > gtol
        iter += 1
        s .= -gvec

        dϕ_0 = dot(s, gvec)
        α, fx = perform_linesearch(ϕ, dϕ, ϕdϕ, 1.0,
                                   fx, dϕ_0, linesearch)
        @. x = x + α*s
        g!(gvec, x)
        gnorm = norm(gvec)
    end

    return (fx, x, iter)
end

# Note that there are many optimization and line search algorithms that allow
# the user to evaluate both the objective and the gradient at the same time, for
# computational efficiency reasons.
# We have included this functionality in the algorithm as the input function `fg!`,
# and even if the Gradient Descent algorithm does not use it explicitly, many of the
# LineSearches algorithms do.

# The Gradient Descent `gdoptimize` method selects a descent direction and calls
# a method `perform_linesearch` that returns the step length `α` and the
# objective value `fx = f(x + α*s)`.
#
# We use multiple dispatch on `linesearch` to call the different line search procedures:

using LineSearches
perform_linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ_0, dϕ_0,
                   linesearch::BackTracking) =
                       linesearch(ϕ, α0, ϕ_0, dϕ_0)
perform_linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ_0, dϕ_0,
                   linesearch::HagerZhang) =
                       linesearch(ϕ, ϕdϕ, α0, ϕ_0, dϕ_0)
perform_linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ_0, dϕ_0,
                   linesearch::MoreThuente) =
                       linesearch(ϕdϕ, α0, ϕ_0, dϕ_0)
perform_linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ_0, dϕ_0,
                   linesearch::StrongWolfe) =
                       linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ_0, dϕ_0)

# The functions ϕ and dϕ represent a univariate objective
# and its derivative, which is used by the line search algorithms.
# To utilize the `fg!` function call in the optimizer, `HagerZhang`,
# `MoreThuente`, and `StrongWolfe` also
# require a function ϕdϕ which returns the univariate objective and the
# derivative at the same time.


# ## Optimizing Rosenbrock
# Here is an example to show how we can combine `gdoptimize` and `LineSearches`
# to minimize the Rosenbrock function, which is defined by

f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

function g!(gvec, x)
    gvec[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
    gvec[2] = 200.0 * (x[2] - x[1]^2)
    gvec
end

function fg!(gvec, x)
    g!(gvec, x)
    f(x)
end


# We can now use `gdoptimize` with `BackTracking` to optimize the Rosenbrock function
# from a given initial condition `x0`.

x0 = [-1., 1.0]
ls = BackTracking(order=3)
fx_bt3, x_bt3, iter_bt3 = gdoptimize(f, g!, fg!, x0, ls)

## Test the results                #src
using Base.Test                    #src
@test fx_bt3 < 1e-12               #src
@test iter_bt3 < 10000             #src
@test x_bt3 ≈ [1.0, 1.0] atol=1e-7 #src

# Interestingly, the `StrongWolfe` line search converges in one iteration, whilst
# all the other algorithms take thousands of iterations.
# This is just luck due to the particular choice of initial condition

ls = StrongWolfe()
fx_sw, x_sw, iter_sw = gdoptimize(f, g!, fg!, x0, ls)

## Test the results               #src
@test fx_sw < 1e-12               #src
@test iter_sw < 10000             #src
@test x_sw ≈ [1.0, 1.0] atol=1e-7 #src


#md # ## [Plain Program](@id customoptimizer-plain-program)
#md #
#md # Below follows a version of the program without any comments.
#md # The file is also available here: [customoptimizer.jl](customoptimizer.jl)
#md #
#md # ```julia
#md # @__CODE__
#md # ```
