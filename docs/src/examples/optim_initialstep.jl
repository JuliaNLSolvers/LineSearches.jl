# # Optim initial step length guess
#
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
using Test                                                #src
@test Optim.f_calls(res_hz) < Optim.f_calls(res_st)       #src
@test Optim.g_calls(res_hz) < Optim.g_calls(res_st)       #src
@test Optim.iterations(res_hz) > Optim.iterations(res_st) #src
