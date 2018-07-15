__precompile__()

module LineSearches

using Printf
import LinearAlgebra: dot, norm
using Parameters, NaNMath

import NLSolversBase
import NLSolversBase: AbstractObjective

export LineSearchException

export BackTracking, HagerZhang, Static, MoreThuente, StrongWolfe

export InitialHagerZhang, InitialStatic, InitialPrevious,
    InitialQuadratic, InitialConstantChange

function make_ϕ(df, x_new, x, s)
    function ϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Evaluate f(x+α*s)
        NLSolversBase.value!(df, x_new)
    end
    ϕ
end
function make_ϕdϕ(df, x_new, x, s)
    function ϕdϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Evaluate ∇f(x+α*s)
        NLSolversBase.value_gradient!(df, x_new)

        # Calculate ϕ(a_i), ϕ'(a_i)
        NLSolversBase.value(df), real(dot(NLSolversBase.gradient(df), s))
    end
    ϕdϕ
end
function make_ϕ_dϕ(df, x_new, x, s)
    function dϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Evaluate ∇f(x+α*s)
        NLSolversBase.gradient!(df, x_new)

        # Calculate ϕ'(a_i)
        real(dot(NLSolversBase.gradient(df), s))
    end
    make_ϕ(df, x_new, x, s), dϕ
end
function make_ϕ_dϕ_ϕdϕ(df, x_new, x, s)
    function dϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Evaluate f(x+α*s) and ∇f(x+α*s)
        NLSolversBase.gradient!(df, x_new)

        # Calculate ϕ'(a_i)
        real(dot(NLSolversBase.gradient(df), s))
    end
    function ϕdϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Evaluate ∇f(x+α*s)
        NLSolversBase.value_gradient!(df, x_new)

        # Calculate ϕ'(a_i)
        NLSolversBase.value(df), real(dot(NLSolversBase.gradient(df), s))
    end
    make_ϕ(df, x_new, x, s), dϕ, ϕdϕ
end
function make_ϕ_ϕdϕ(df, x_new, x, s)
    function ϕdϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Evaluate ∇f(x+α*s)
        NLSolversBase.value_gradient!(df, x_new)

        # Calculate ϕ'(a_i)
        NLSolversBase.value(df), real(dot(NLSolversBase.gradient(df), s))
    end
    make_ϕ(df, x_new, x, s), ϕdϕ
end

include("types.jl")

# Line Search Methods
include("backtracking.jl")
include("strongwolfe.jl")
include("morethuente.jl")
include("hagerzhang.jl") # Also includes InitialHagerZhang
include("static.jl")

# Initial guess methods
include("initialguess.jl")

include("deprecate.jl")

end # module
