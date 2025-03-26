# https://github.com/JuliaNLSolvers/LineSearches.jl/issues/151
using LineSearches, LinearAlgebra, Test

A = randn(100, 100)
x0 = randn(100)
b = A * x0

# Objective function and gradient
f(x) = 0.5 * norm(A * x - b)^2
g!(gvec, x) = (gvec .= A' * (A * x - b))
fg!(gvec, x) = (g!(gvec, x); return f(x))

# Init
x = 1.0f1 * randn(100)
gv = similar(x)

# Line search
α0 = 1.0f-3
ϕ0 = fg!(gv, x)
s = -1 * gv
dϕ0 = dot(gv, s)
println(ϕ0, ", ", dϕ0)

# Univariate line search functions
ϕ(α) = f(x .+ α .* s)
function dϕ(α)
    g!(gv, x .+ α .* s)
    return dot(gv, s)
end
function ϕdϕ(α)
    phi = fg!(gv, x .+ α .* s)
    dphi = dot(gv, s)
    return (phi, dphi)
end

res = (StrongWolfe())(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
@test res[2] > 0
@test res[2] == ϕ(res[1])

struct LineSearchTestCase
    alphas::Vector{Float64}
    values::Vector{Float64}
    slopes::Vector{Float64}
end

@testset "HZ convergence issues" begin
    @testset "Flatness check issues" begin
        function prepare_test_case(; alphas, values, slopes)
            perm = sortperm(alphas)
            alphas = alphas[perm]
            push!(alphas, alphas[end] + 1)
            values = values[perm]
            push!(values, values[end])
            slopes = slopes[perm]
            push!(slopes, 0.0)
            return LineSearchTestCase(alphas, values, slopes)
        end

        tc1 = prepare_test_case(;
            alphas = [0.0, 1.0, 5.0, 3.541670844449739],
            values = [
                3003.592409634743,
                2962.0378569864743,
                2891.4462095232184,
                3000.9760725116876,
            ],
            slopes = [
                -22332.321416890798,
                -20423.214551925797,
                11718.185026267562,
                -22286.821227217057,
            ],
        )

        function tc_to_f(tc)
            function f(x)
                i = findfirst(u -> u > x, tc.alphas) - 1
                xk = tc.alphas[i]
                xkp1 = tc.alphas[i+1]
                dx = xkp1 - xk
                t = (x - xk) / dx
                h00t = 2t^3 - 3t^2 + 1
                h10t = t * (1 - t)^2
                h01t = t^2 * (3 - 2t)
                h11t = t^2 * (t - 1)
                val =
                    h00t * tc.values[i] +
                    h10t * dx * tc.slopes[i] +
                    h01t * tc.values[i+1] +
                    h11t * dx * tc.slopes[i+1]

                return val
            end
        end
        function tc_to_fdf(tc)
            function fdf(x)
                i = findfirst(u -> u > x, tc.alphas) - 1
                xk = tc.alphas[i]
                xkp1 = tc.alphas[i+1]
                dx = xkp1 - xk
                t = (x - xk) / dx
                h00t = 2t^3 - 3t^2 + 1
                h10t = t * (1 - t)^2
                h01t = t^2 * (3 - 2t)
                h11t = t^2 * (t - 1)
                val =
                    h00t * tc.values[i] +
                    h10t * dx * tc.slopes[i] +
                    h01t * tc.values[i+1] +
                    h11t * dx * tc.slopes[i+1]

                h00tp = 6t^2 - 6t
                h10tp = 3t^2 - 4t + 1
                h01tp = -6t^2 + 6 * t
                h11tp = 3t^2 - 2t
                slope =
                    (
                        h00tp * tc.values[i] +
                        h10tp * dx * tc.slopes[i] +
                        h01tp * tc.values[i+1] +
                        h11tp * dx * tc.slopes[i+1]
                    ) / dx
                println(x, " ", val, " ", slope)
                return val, slope
            end
        end

        function test_tc(tc, check_flatness)
            cache = LineSearchCache{Float64}()
            hz = HagerZhang(; cache, check_flatness)
            f = tc_to_f(tc)
            fdf = tc_to_fdf(tc)
            hz(f, fdf, 1.0, fdf(0.0)...), cache
        end

        res, res_cache = test_tc(tc1, true)
        @show res
        @show res_cache
        @test_broken minimum(res_cache.values) == res[2]

        res2, res_cache2 = test_tc(tc1, false)
        @test minimum(res_cache2.values) == res2[2]
        #=
        using AlgebraOfGraphics, CairoMakie
        draw(data((x=0.0:0.05:5.5, y=map(x->tc_to_f(tc1)(x), 0:0.05:5.5)))*mapping(:x,:y)*visual(Scatter)+
        data((alphas=res_cache.alphas, values=res_cache.values))*mapping(:alphas,:values)*visual(Scatter; color=:red))
        =#
    end

    # should add as upstream
    #=
    @testset "from kbarros" begin
        # The minimizer is x0=[0, 2πn/100], with f(x0) = 1. Any integer n is fine.
        function f(x)
            return (x[1]^2 + 1) * (2 - cos(100*x[2]))
        end

        using Optim

        function test_converges(method)
            for i in 1:100
                r = randn(2)
                res = optimize(f, r, method)
                if Optim.converged(res) && minimum(res) > f([0,0]) + 1e-8
                    println("""
                        Incorrectly reported convergence after $(res.iterations) iterations
                        Reached x = $(Optim.minimizer(res)) with f(x) = $(minimum(res))
                        """)
                end
            end
        end

        # Works successfully, no printed output
        test_converges(LBFGS(; linesearch=Optim.LineSearches.BackTracking(order=2)))

        # Prints ~10 failures to converge (in 100 tries). Frequently fails after the
        # first line search.
        test_converges(ConjugateGradient(; linesearch=Optim.LineSearches.HagerZhang(check_flatness=false)))
    end
    =#
end
