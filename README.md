# LineSearches

[![Build Status](https://travis-ci.org/JuliaNLSolvers/LineSearches.jl.svg?branch=master)](https://travis-ci.org/JuliaNLSolvers/LineSearches.jl)
[![Codecov branch](https://img.shields.io/codecov/c/github/JuliaNLSolvers/LineSearches.jl/master.svg?maxAge=2592000)](https://codecov.io/gh/JuliaNLSolvers/LineSearches.jl)

## Description
This package provides an interface to line search algorithms implemented in Julia.
The code was originally written as part of [Optim](https://github.com/JuliaNLSolvers/Optim.jl),
but has now been separated out to its own package.

### Available line search algorithms
In Example 1 we show how to choose between the line search algorithms
in `Optim`.
* `HagerZhang` (Taken from the Conjugate Gradient implementation
  by Hager and Zhang, 2006)
* `MoreThuente` (From the algorithm in More and Thuente, 1994)
* `BackTracking` (Described in Nocedal and Wright, 2006)
* `StrongWolfe` (Nocedal and Wright)
* `Static` (Simply takes a given step length, default 1.0)

### Available initial step length procedures
The package provides some procedures to calculate the initial step
length that is passed to the line search algorithm. See Example 2 for
its usage in `Optim`.
* `InitialPrevious` (Use the step length from the previous
  optimization iteration)
* `InitialStatic` (Use the same initial step length each time)
* `InitialHagerZhang` (Taken from Hager and Zhang, 2006)
* `InitialQuadratic` (Propose initial step length based on a quadratic
  interpolation)
* `InitialConstantChange` (Propose initial step length assuming
  constant change in step length)


## Example 1
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
   * |x - x'| ≤ 1.0e-32: false 
     |x - x'| = 3.06e-09 
   * |f(x) - f(x')| ≤ 1.0e-32 |f(x)|: false
     |f(x) - f(x')| = 3.03e+13 |f(x)|
   * |g(x)| ≤ 1.0e-08: true 
     |g(x)| = 1.11e-15 
   * Stopped by an increasing objective: false
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
   * |x - x'| ≤ 1.0e-32: false 
     |x - x'| = 1.36e-05 
   * |f(x) - f(x')| ≤ 1.0e-32 |f(x)|: false
     |f(x) - f(x')| = 1.21e+08 |f(x)|
   * |g(x)| ≤ 1.0e-08: true 
     |g(x)| = 4.16e-09 
   * Stopped by an increasing objective: false
   * Reached Maximum Number of Iterations: false
 * Objective Calls: 19
 * Gradient Calls: 15
 * Hessian Calls: 14
```

## Example 2
This example shows how to use the initial step length procedures with `Optim`.
We solve the Rosenbrock problem with two different procedures.

First, run `Newton` with the (default) initial guess and line search procedures.
```julia
using Optim, LineSearches
prob = Optim.UnconstrainedProblems.examples["Rosenbrock"]

algo_st = Newton(alphaguess = InitialStatic(), linesearch = HagerZhang())
res_st = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_st)
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
   * |x - x'| ≤ 1.0e-32: false 
     |x - x'| = 3.06e-09 
   * |f(x) - f(x')| ≤ 1.0e-32 |f(x)|: false
     |f(x) - f(x')| = 3.03e+13 |f(x)|
   * |g(x)| ≤ 1.0e-08: true 
     |g(x)| = 1.11e-15 
   * Stopped by an increasing objective: false
   * Reached Maximum Number of Iterations: false
 * Objective Calls: 44
 * Gradient Calls: 44
 * Hessian Calls: 14
```

We can now try with the initial step length guess from Hager and Zhang.
``` julia
algo_prev = Newton(alphaguess = InitialHagerZhang(α0=1.0), linesearch = HagerZhang())
res_prev = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_prev)
```

This gives the following result, reducing the number of function and gradient calls, but increasing the number of iterations.
``` julia
Results of Optimization Algorithm
 * Algorithm: Newton's Method
 * Starting Point: [0.0,0.0]
 * Minimizer: [0.9999999974436653,0.9999999948855858]
 * Minimum: 6.535152e-18
 * Iterations: 15
 * Convergence: true
   * |x - x'| ≤ 1.0e-32: false 
     |x - x'| = 1.09e-05 
   * |f(x) - f(x')| ≤ 1.0e-32 |f(x)|: false
     |f(x) - f(x')| = 8.61e+08 |f(x)|
   * |g(x)| ≤ 1.0e-08: true 
     |g(x)| = 4.41e-09 
   * Stopped by an increasing objective: false
   * Reached Maximum Number of Iterations: false
 * Objective Calls: 36
 * Gradient Calls: 21
 * Hessian Calls: 15
```


## References
- W. W. Hager and H. Zhang (2006) "Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed descent." ACM Transactions on Mathematical Software 32: 113-137.
- Moré, Jorge J., and David J. Thuente. "Line search algorithms with guaranteed sufficient decrease." ACM Transactions on Mathematical Software (TOMS) 20.3 (1994): 286-307.
- Nocedal, Jorge, and Stephen Wright. "Numerical optimization." Springer Science & Business Media, 2006.
