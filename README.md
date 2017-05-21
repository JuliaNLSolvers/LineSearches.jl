# LineSearches

[![Build Status](https://travis-ci.org/JuliaNLSolvers/LineSearches.jl.svg?branch=master)](https://travis-ci.org/JuliaNLSolvers/LineSearches.jl)
[![Codecov branch](https://img.shields.io/codecov/c/github/JuliaNLSolvers/LineSearches.jl/master.svg?maxAge=2592000)](https://codecov.io/gh/JuliaNLSolvers/LineSearches.jl)

## Description
This package provides an interface to line search algorithms implemented in Julia.
The code was originally written as part of [Optim](https://github.com/JuliaNLSolvers/Optim.jl),
but has now been separated out to its own package.

### Available line search algorithms
* `HagerZhang` (Taken from the Conjugate Gradient implementation
  by Hager and Zhang, 2006)
* `MoreThuente` (From the algorithm in More and Thuente, 1994)
* `BackTracking` (Described in Nocedal and Wright, 2006)
* `StrongWolfe` (Nocedal and Wright)
* `Static` (Simply takes a given step length, default 1.0)

## Example
This example shows how to use `LineSearches` with `Optim`.
We solve the Rosenbrock problem with two different line search algorithms.

First, run `Newton` with the default line search algorithm:
```julia
using Optim, LineSearches
prob = Optim.UnconstrainedProblems.examples["Rosenbrock"]

algo_hz = Newton(linesearch = HagerZhang())
res_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)
```

This gives the result
``` julia
Results of Optimization Algorithm
 * Algorithm: Newton's Method
 * Starting Point: [0.0,0.0]
 * Minimizer: [0.9999999999999994,0.9999999999999989]
 * Minimum: 3.081488e-31
 * Iterations: 14
 * Convergence: true
   * |x - x'| < 1.0e-32: false
   * |f(x) - f(x')| / |f(x)| < 1.0e-32: false
   * |g(x)| < 1.0e-08: true
   * f(x) > f(x'): false
   * Reached Maximum Number of Iterations: false
 * Objective Calls: 44
 * Gradient Calls: 44
 * Hessian Calls: 14
```

Now we can try `Newton` with the cubic backtracking line search:
``` julia
algo_bt3 = Newton(linesearch = BackTracking(order=3))
res_bt3 = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_bt3)
```

This gives the following result, reducing the number of function and gradient calls:
``` julia
Results of Optimization Algorithm
 * Algorithm: Newton's Method
 * Starting Point: [0.0,0.0]
 * Minimizer: [0.9999999959215587,0.9999999918223065]
 * Minimum: 1.667699e-17
 * Iterations: 14
 * Convergence: true
   * |x - x'| < 1.0e-32: false
   * |f(x) - f(x')| / |f(x)| < 1.0e-32: false
   * |g(x)| < 1.0e-08: true
   * f(x) > f(x'): false
   * Reached Maximum Number of Iterations: false
 * Objective Calls: 19
 * Gradient Calls: 15
 * Hessian Calls: 14
```

## References
- W. W. Hager and H. Zhang (2006) "Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed descent." ACM Transactions on Mathematical Software 32: 113-137.
- Moré, Jorge J., and David J. Thuente. "Line search algorithms with guaranteed sufficient decrease." ACM Transactions on Mathematical Software (TOMS) 20.3 (1994): 286-307.
- Nocedal, Jorge, and Stephen Wright. "Numerical optimization." Springer Science & Business Media, 2006.
