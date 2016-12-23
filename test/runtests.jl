using LineSearches
using Base.Test

lsfunctions = (hagerzhang!, strongwolfe!,
               morethuente!, backtracking!,
               bt2!, bt3!)

println("Running tests:")
my_tests = [
    "alphacalc.jl",
    "optim.jl",
    "backtracking.jl"
]

for my_test in my_tests
    println("\n * $(my_test)")
    include(my_test)
end
