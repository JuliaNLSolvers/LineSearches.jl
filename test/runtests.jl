using LineSearches
using   Compat,
        Compat.Test,
        Compat.LinearAlgebra
using OptimTestProblems
if Sys.WORD_SIZE != 32
    # Bug in HigherPrecision, waiting for
    # https://github.com/saschatimme/HigherPrecision.jl/pull/21
    using HigherPrecision
end
import NLSolversBase

debug_printing = false

lstypes =  (Static(), HagerZhang(), StrongWolfe(), MoreThuente(),
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
