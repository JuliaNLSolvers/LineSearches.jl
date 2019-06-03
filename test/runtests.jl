using LineSearches
using Test
using LinearAlgebra: norm, dot
using OptimTestProblems
using DoubleFloats: Double64, Double32, Double16
import NLSolversBase

debug_printing = false

lstypes =  (Static(), HagerZhang(), StrongWolfe(), MoreThuente(),
            BackTracking(), BackTracking(order=2) )

my_tests = [
    "initial.jl",
    "alphacalc.jl",
    "arbitrary_precision.jl",
    "examples.jl"
]

mutable struct StateDummy
    alpha
    x
    x_ls
    f_x_previous
    s
end

for my_test in my_tests
    println("\n * $(my_test)")
    @time include(my_test)
end

@testset "External constructors" begin
    bt = BackTracking{Float32}()
    @test typeof(bt) == BackTracking{Float32, Int}

    hz = HagerZhang{Float32}()
    @test typeof(hz) == HagerZhang{Float32, Base.RefValue{Bool}}
end
