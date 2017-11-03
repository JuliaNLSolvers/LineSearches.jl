using LineSearches
using Base.Test
using OptimTestProblems

debug_printing = false

lstypes =  (Static(), Static(scaled=true), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

my_tests = [
    "initial.jl",
    "api.jl",
    "alphacalc.jl"
]

mutable struct StateDummy
    alpha
    x
    x_ls
    f_x_previous
    s
    lsr::LineSearchResults
    mayterminate::Bool
    dphi0_previous
end

for my_test in my_tests
    println("\n * $(my_test)")
    @time include(my_test)
end
