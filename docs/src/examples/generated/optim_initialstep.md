```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/LineSearches/docs/src/examples/optim_initialstep.jl"
```

# Optim initial step length guess

!!! tip
    This example is also available as a Jupyter notebook:
    [`optim_initialstep.ipynb`](https://nbviewer.jupyter.org/github/TRAVIS_REPO_SLUG/blob/gh-pages/TRAVIS_TAG/examples/generated/optim_initialstep.ipynb)

This example shows how to use the initial step length procedures
with [Optim](https://github.com/JuliaNLSolvers/Optim.jl).  We solve
the Rosenbrock problem with two different procedures.

First, run `Newton` with the (default) initial guess and line search procedures.

```@example optim_initialstep
using Optim, LineSearches
import OptimTestProblems.MultivariateProblems
UP = MultivariateProblems.UnconstrainedProblems
prob = UP.examples["Rosenbrock"]

algo_st = Newton(alphaguess = InitialStatic(), linesearch = HagerZhang())
res_st = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_st)
```

We can now try with the initial step length guess from Hager and Zhang.

```@example optim_initialstep
algo_hz = Newton(alphaguess = InitialHagerZhang(α0=1.0), linesearch = HagerZhang())
res_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)
```

From the result we see that this has reduced the number of function and gradient calls, but increased the number of iterations.

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

