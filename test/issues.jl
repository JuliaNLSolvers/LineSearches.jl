# https://github.com/JuliaNLSolvers/LineSearches.jl/issues/151
using LineSearches, LinearAlgebra, Test

A = randn(100, 100)
x0 = randn(100)
b = A*x0

# Objective function and gradient
f(x) = .5*norm(A*x - b)^2
g!(gvec, x) = (gvec .= A'*(A*x-b))
fg!(gvec, x) = (g!(gvec, x); return f(x))

# Init
x = 1f1*randn(100)
gv = similar(x)

# Line search
α0 = 1f-3
ϕ0 = fg!(gv, x)
s = -1*gv
dϕ0 = dot(gv, s)
# println(ϕ0, ", ", dϕ0)

# Univariate line search functions
ϕ(α) = f(x .+ α.*s)
function dϕ(α)
    g!(gv, x .+ α.*s)
    return dot(gv, s)
end
function ϕdϕ(α)
    phi = fg!(gv, x .+ α.*s)
    dphi = dot(gv, s)
    return (phi, dphi)
end

res = (StrongWolfe())(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
@test res[2] > 0
@test res[2] == ϕ(res[1])

# Flatness check in HagerZhang
function makeϕdϕ(a)
    @assert axes(a) == (1:2,)
    A = a*a'
    f(x) = x'*A*x/2   # Hessian has one positive and one zero eigenvalue (in exact arithmetic)
    df(x) = A*x
    x0 = [a[2], -a[1]]
    d = -x0 / 2
    ϕ(α) = f(x0 + α*d)
    dϕ(α) = dot(df(x0 + α*d), d)
    ϕdϕ(α) = (ϕ(α), dϕ(α))
    return ϕ, dϕ, ϕdϕ
end

@testset "Flatness" begin
    cache = LineSearchCache{Float64}()
    lsalgs =  (HagerZhang(; cache=cache), StrongWolfe(; cache=cache), MoreThuente(; cache=cache),
               BackTracking(; cache=cache), BackTracking(; order=2, cache=cache) )

    npass = zeros(Int, length(lsalgs))
    n, nmax = 0, 1000
    while n < nmax
        ϕ, dϕ, ϕdϕ = makeϕdϕ(randn(2))
        ϕ0, dϕ0 = ϕdϕ(0)
        dϕ0 < -eps(abs(ϕ0)) || continue    # any "slope" is just roundoff error, but we want roundoff that looks like descent
        n += 1
        for (i, ls) in enumerate(lsalgs)
            res = ls(ϕ, dϕ, ϕdϕ, 1.0, ϕ(0.0), dϕ(0.0))
            npass[i] += length(cache.alphas) < 10
        end
    end
    @test_broken all(npass .== nmax)
    @show npass
    @test all(npass[[3, 4, 5]] .>= nmax-3)
end
