using Optim, LineSearches
import OptimTestProblems.MultivariateProblems
UP = MultivariateProblems.UnconstrainedProblems
prob = UP.examples["Rosenbrock"]

algo_st = Newton(alphaguess = InitialStatic(), linesearch = HagerZhang())
res_st = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_st)

algo_hz = Newton(alphaguess = InitialHagerZhang(α0=1.0), linesearch = HagerZhang())
res_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

