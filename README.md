# LineSearches

[![Build Status](https://travis-ci.org/JuliaNLSolvers/LineSearches.jl.svg?branch=master)](https://travis-ci.org/JuliaNLSolvers/LineSearches.jl)
[![Codecov branch](https://img.shields.io/codecov/c/github/JuliaNLSolvers/LineSearches.jl/master.svg?maxAge=2592000)](https://codecov.io/gh/JuliaNLSolvers/LineSearches.jl)
[![][docs-stable-img]][docs-stable-url]

## Description
This package provides an interface to line search algorithms implemented in Julia.
The code was originally written as part of [Optim](https://github.com/JuliaNLSolvers/Optim.jl),
but has now been separated out to its own package.

### Available line search algorithms
In [the docs](https://julianlsolvers.github.io/LineSearches.jl/latest/examples/generated/optim_linesearch.html) we show how to choose between the line search algorithms
in `Optim`.
* `HagerZhang` (Taken from the Conjugate Gradient implementation
  by Hager and Zhang, 2006)
* `MoreThuente` (From the algorithm in More and Thuente, 1994)
* `BackTracking` (Described in Nocedal and Wright, 2006)
* `StrongWolfe` (Nocedal and Wright)
* `Static` (Takes the proposed initial step length.)

### Available initial step length procedures
The package provides some procedures to calculate the initial step
length that is passed to the line search algorithm. See [the docs](https://julianlsolvers.github.io/LineSearches.jl/latest/examples/generated/optim_initialstep.html) for
its usage in `Optim`.
* `InitialPrevious` (Use the step length from the previous
  optimization iteration)
* `InitialStatic` (Use the same initial step length each time)
* `InitialHagerZhang` (Taken from Hager and Zhang, 2006)
* `InitialQuadratic` (Propose initial step length based on a quadratic
  interpolation)
* `InitialConstantChange` (Propose initial step length assuming
  constant change in step length)


## Documentation
For more details and options, see the documentation
- [STABLE][docs-stable-url] — most recently tagged version of the documentation.
- [LATEST][docs-latest-url] — in-development version of the documentation.

## Example usage
Here is how to get a simple linesearch for a one-dimensional function working:
```julia
using LineSearches

ϕ(x) = (x - π)^4
dϕ(x) = 4*(x-π)^3
ϕdϕ(x) = ϕ(x),dϕ(x)

α0 = 9.0
ϕ0 = ϕ(0.0)
dϕ0 = dϕ(0.0)

for ls in (Static,BackTracking,HagerZhang,MoreThuente,StrongWolfe)
    res = (ls())(ϕ, dϕ, ϕdϕ, α0, ϕ0,dϕ0)
    println(ls, ": ", res)
end
```
For more examples, see the documentation.

## References
- W. W. Hager and H. Zhang (2006) "Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed descent." ACM Transactions on Mathematical Software 32: 113-137.
- Moré, Jorge J., and David J. Thuente. "Line search algorithms with guaranteed sufficient decrease." ACM Transactions on Mathematical Software (TOMS) 20.3 (1994): 286-307.
- Nocedal, Jorge, and Stephen Wright. "Numerical optimization." Springer Science & Business Media, 2006.


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://julianlsolvers.github.io/LineSearches.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://julianlsolvers.github.io/LineSearches.jl/stable
