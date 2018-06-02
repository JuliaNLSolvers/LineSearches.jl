# # Optim initial step length guess
#
#src TODO: Find a way to run these with Literate when deploying via Travis
#src TODO: This file must currently be run locally and not on CI, and then
#src TODO: the md file must be copied over to the correct directory.
#src TODO: The reason is that there may be breaking changes between Optim and LineSearches,
#src TODO: so we don't want that to mess up JuliaCIBot
#-
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`optim_initialstep.ipynb`](@__NBVIEWER_ROOT_URL__examples/generated/optim_initialstep.ipynb)
#-
#
# This example shows how to use the initial step length procedures
# with [Optim](https://github.com/JuliaNLSolvers/Optim.jl).  We solve
# the Rosenbrock problem with two different procedures.
#
# First, run `Newton` with the (default) initial guess and line search procedures.
using Optim, LineSearches
import OptimTestProblems.MultivariateProblems
UP = MultivariateProblems.UnconstrainedProblems
prob = UP.examples["Rosenbrock"]

algo_st = Newton(alphaguess = InitialStatic(), linesearch = HagerZhang())
res_st = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_st)


# We can now try with the initial step length guess from Hager and Zhang.
algo_hz = Newton(alphaguess = InitialHagerZhang(α0=1.0), linesearch = HagerZhang())
res_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)


# From the result we see that this has reduced the number of function and gradient calls, but increased the number of iterations.

## Test the results                                       #src
using Base.Test                                           #src
@test Optim.f_calls(res_hz) < Optim.f_calls(res_st)       #src
@test Optim.g_calls(res_hz) < Optim.g_calls(res_st)       #src
@test Optim.iterations(res_hz) > Optim.iterations(res_st) #src
