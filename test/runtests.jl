using LineSearches
using Base.Test

debug_printing = true

lstypes =  (Static(), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

dep_lsfunctions = (basic!, hagerzhang!, strongwolfe!, morethuente!,
                   backtracking!, bt3!, bt2!)

my_tests = [
    "alphacalc.jl",
    "optim.jl",
    "backtracking.jl",
    "counter.jl"
]

for my_test in my_tests
    println("\n * $(my_test)")
    @time include(my_test)
end
