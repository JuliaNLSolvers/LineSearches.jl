# LineSearches

[![Build Status](https://travis-ci.org/anriseth/LineSearches.jl.svg?branch=master)](https://travis-ci.org/anriseth/LineSearches.jl)
[![Codecov branch](https://img.shields.io/codecov/c/github/anriseth/LineSearches.jl/master.svg?maxAge=2592000)](https://codecov.io/gh/anriseth/LineSearches.jl)

## Description
This package provides an interface to line search algorithms implemented in Julia.
The code was originally written as part of [Optim](https://github.com/JuliaOpt/Optim.jl),
but has now been separated out to its own package.

### Available line search algorithms
* `hagerzhang!` (Taken from the Conjugate Gradient implementation
  by Hager and Zhang (HZ))
* `backtracking!`
* `interpbacktracking!`
* `strongwolfe!`
* `morethuente!`


## Example
This example shows how to use `LineSearches` with `Optim`.
We solve the Rosenbrock problem with two different line search algorithms.

First, run `Newton` with the default line search algorithm:
```julia
using Optim
prob = Optim.UnconstrainedProblems.examples["Rosenbrock"]

algo_hz = Newton(;linesearch! = hagerzhang!)
res_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)
```

This gives the result
``` julia
Results of Optimization Algorithm
 * Algorithm: Newton's Method
 * Starting Point: [0.0,0.0]
 * Minimizer: [0.9999999999979515,0.9999999999960232]
 * Minimum: 5.639268e-24
 * Iterations: 13
 * Convergence: true
   * |x - x'| < 1.0e-32: false
   * |f(x) - f(x')| / |f(x)| < 1.0e-32: false
   * |g(x)| < 1.0e-08: true
   * Reached Maximum Number of Iterations: false
 * Objective Function Calls: 54
 * Gradient Calls: 54
```

Now we can try `Newton` with the Mure Thuente line search:
``` julia
algo_mt = Newton(;linesearch! = morethuente!)
res_mt = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_mt)
```

This gives the following result, reducing the number of function and gradient calls:
``` julia
 * Algorithm: Newton's Method
 * Starting Point: [0.0,0.0]
 * Minimizer: [0.9999999999999992,0.999999999999998]
 * Minimum: 2.032549e-29
 * Iterations: 14
 * Convergence: true
   * |x - x'| < 1.0e-32: false
   * |f(x) - f(x')| / |f(x)| < 1.0e-32: false
   * |g(x)| < 1.0e-08: true
   * Reached Maximum Number of Iterations: false
 * Objective Function Calls: 45
 * Gradient Calls: 45
```

## References
W. W. Hager and H. Zhang (2006) "Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed descent." ACM Transactions on Mathematical Software 32: 113-137.
