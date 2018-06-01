```@meta
DocTestSetup = :(using LineSearches)
```

# LineSearches.jl
*A line search toolbox written in Julia.*

## Introduction
`LineSearches` provides a collection of line search routines for
optimization and nonlinear solvers.
The package can be used on its own, but it also provides extra
supporting functionality for `Optim.jl` and `NLsolve.jl`.


## Installation

To install, simply run the following in the Julia REPL:
```julia
Pkg.add("LineSearches")
```
and then run
```julia
using JuAFEM
```
to load the package.
