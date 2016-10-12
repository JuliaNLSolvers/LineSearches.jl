using LineSearches
using Base.Test
import Optim

lsfunctions = (hz_linesearch!, interpolating_linesearch!,
               mt_linesearch!, backtracking_linesearch!,
               interpbacktrack_linesearch!)

println("Running tests:")
my_tests = [
    "alphacalc.jl",
#    "optim.jl"
]

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end
