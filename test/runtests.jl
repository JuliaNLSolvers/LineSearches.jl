using LineSearches
using Base.Test
import Optim

lsfunctions = (hagerzhang!, strongwolfe!,
               morethuente!, backtracking!,
               interpbacktrack!)

println("Running tests:")
my_tests = [
    "alphacalc.jl",
#    "optim.jl"
]

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end
