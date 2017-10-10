using LineSearches
using Base.Test

debug_printing = false

lstypes =  (Static(), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

my_tests = [
    "alphacalc.jl",
    "backtracking.jl"
]

for my_test in my_tests
    println("\n * $(my_test)")
    @time include(my_test)
end
