using LineSearches
using   Compat,
        Test,
        Compat.LinearAlgebra
using OptimTestProblems
using HigherPrecision
import NLSolversBase

debug_printing = false

lstypes =  (Static(), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

my_tests = [
    "initial.jl",
    "alphacalc.jl",
    "arbitrary_precision.jl",
    "examples.jl"
]

mutable struct StateDummy
    alpha
    x
    x_ls
    f_x_previous
    s
end

for my_test in my_tests
    println("\n * $(my_test)")
    @time include(my_test)
end

# Build the docs
if get(ENV, "TRAVIS_OS_NAME", "") == "linux" && get(ENV, "TRAVIS_JULIA_VERSION", "") == "0.6"
    include(joinpath(@__DIR__, "../docs/make.jl"))
end
