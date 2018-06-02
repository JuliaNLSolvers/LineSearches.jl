```@meta
DocTestSetup = :(using LineSearches)
```

# LineSearches.jl
*A line search toolbox written in Julia.*

## Introduction
`LineSearches` provides a collection of line search routines for
optimization and nonlinear solvers.  The package can be used on its
own, but it also provides extra supporting functionality for
[Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) and
[NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl).


## Available line search algorithms
* `HagerZhang` (Taken from the Conjugate Gradient implementation
  by Hager and Zhang, 2006)
* `MoreThuente` (From the algorithm in More and Thuente, 1994)
* `BackTracking` (Described in Nocedal and Wright, 2006)
* `StrongWolfe` (Nocedal and Wright)
* `Static` (Takes the proposed initial step length.)

## Available initial step length procedures
The package provides some procedures to calculate the initial step
length that is passed to the line search algorithm, currently specialized to
be used with Optim and NLsolve.
* `InitialPrevious` (Use the step length from the previous
  optimization iteration)
* `InitialStatic` (Use the same initial step length each time)
* `InitialHagerZhang` (Taken from Hager and Zhang, 2006)
* `InitialQuadratic` (Propose initial step length based on a quadratic
  interpolation)
* `InitialConstantChange` (Propose initial step length assuming
  constant change in step length)

## Installation

To install, simply run the following in the Julia REPL:
```julia
Pkg.add("LineSearches")
```
and then run
```julia
using LineSearches
```
to load the package.


## References

- W. W. Hager and H. Zhang (2006) "Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed descent." ACM Transactions on Mathematical Software 32: 113-137.
- Moré, Jorge J., and David J. Thuente. "Line search algorithms with guaranteed sufficient decrease." ACM Transactions on Mathematical Software (TOMS) 20.3 (1994): 286-307.
- Nocedal, Jorge, and Stephen Wright. "Numerical optimization." Springer Science & Business Media, 2006.
