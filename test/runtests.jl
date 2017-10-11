using LineSearches
using Base.Test
using OptimTestProblems

debug_printing = false

lstypes =  (Static(), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

my_tests = [
    "api.jl",
    "alphacalc.jl",
    "backtracking.jl"
]

for my_test in my_tests
    println("\n * $(my_test)")
    @time include(my_test)
end
