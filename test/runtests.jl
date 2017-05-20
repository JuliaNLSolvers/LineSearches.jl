using LineSearches
using Base.Test

lstypes =  (Static(), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

dep_lsfunctions = (basic!, hagerzhang!, strongwolfe!, morethuente!,
                  backtracking!, bt2!, bt3!)

println("Running tests:")
my_tests = [
    "alphacalc.jl",
# TODO: Add Optim tests back when Optim.jl#388 is merged
#    "optim.jl",
#    "backtracking.jl",
#    "counter.jl"
]

for my_test in my_tests
    println("\n * $(my_test)")
    include(my_test)
end
