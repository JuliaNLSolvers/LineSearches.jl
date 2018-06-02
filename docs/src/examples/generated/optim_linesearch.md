```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/LineSearches/docs/src/examples/optim_linesearch.jl"
```

# Optim line search

!!! tip
    This example is also available as a Jupyter notebook:
    [`optim_linesearch.ipynb`](https://nbviewer.jupyter.org/github/TRAVIS_REPO_SLUG/blob/gh-pages/TRAVIS_TAG/examples/generated/optim_linesearch.ipynb)

This example shows how to use `LineSearches` with
[Optim](https://github.com/JuliaNLSolvers/Optim.jl).  We solve the
Rosenbrock problem with two different line search algorithms.

First, run `Newton` with the default line search algorithm:

```@example optim_linesearch
using Optim, LineSearches
import OptimTestProblems.MultivariateProblems
UP = MultivariateProblems.UnconstrainedProblems
prob = UP.examples["Rosenbrock"]

algo_hz = Newton(linesearch = HagerZhang())
res_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)
```

Now we can try `Newton` with the cubic backtracking line search,
which reduced the number of objective and gradient calls.

```@example optim_linesearch
algo_bt3 = Newton(linesearch = BackTracking(order=3))
res_bt3 = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_bt3)
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

