using LineSearches
using   Compat,
        Compat.Test,
        Compat.LinearAlgebra
using OptimTestProblems
using HigherPrecision
import NLSolversBase

debug_printing = false

lstypes =  (Static(), Static(scaled=true), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

my_tests = [
    "initial.jl",
    "alphacalc.jl",
    "arbitrary_precision.jl"
]

mutable struct StateDummy
    alpha
    x
    x_ls
    f_x_previous
    s
    dphi_0_previous
end

for my_test in my_tests
    println("\n * $(my_test)")
    @time include(my_test)
end
