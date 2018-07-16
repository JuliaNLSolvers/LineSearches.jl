# # Optim line search
#
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`optim_linesearch.ipynb`](@__NBVIEWER_ROOT_URL__examples/generated/optim_linesearch.ipynb)
#-
#
# This example shows how to use `LineSearches` with
# [Optim](https://github.com/JuliaNLSolvers/Optim.jl).  We solve the
# Rosenbrock problem with two different line search algorithms.
#
# First, run `Newton` with the default line search algorithm:

using Optim, LineSearches
import OptimTestProblems.MultivariateProblems
UP = MultivariateProblems.UnconstrainedProblems
prob = UP.examples["Rosenbrock"]

algo_hz = Newton(linesearch = HagerZhang())
res_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)

# Now we can try `Newton` with the cubic backtracking line search,
# which reduced the number of objective and gradient calls.

algo_bt3 = Newton(linesearch = BackTracking(order=3))
res_bt3 = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_bt3)


## Test the results                                  #src
using Test                                           #src
@test Optim.f_calls(res_bt3) < Optim.f_calls(res_hz) #src
@test Optim.g_calls(res_bt3) < Optim.g_calls(res_hz) #src
