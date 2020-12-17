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
println(ϕ0, ", ", dϕ0)

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
