module TestCases

# A module for reproducing simplified versions of objective functions that cause trouble for line search methods

using NLSolversBase

export LineSearchTestCase

struct LineSearchTestCase
    alphas::Vector{Float64}
    values::Vector{Float64}
    slopes::Vector{Float64}

    function LineSearchTestCase(alphas, values, slopes)
        Base.require_one_based_indexing(alphas, values, slopes)
        n = length(alphas)
        if n != length(values) || n != length(slopes)
            throw(ArgumentError("Lengths of alphas, values, and slopes must match"))
        end
        # Ensure ordered & unique
        perm = sortperm(alphas)
        delidx = Int[]
        for i in firstindex(perm)+1:lastindex(perm)
            if alphas[perm[i]] == alphas[perm[i-1]]
                push!(delidx, i)
            end
        end
        deleteat!(perm, delidx)
        alphas, values, slopes = alphas[perm], values[perm], slopes[perm]
        # For interpolation, add a dummy value at the end
        push!(alphas, alphas[end]+1)
        push!(values, values[end])
        push!(slopes, 0.0)
        return new(alphas, values, slopes)
    end
end

Base.minimum(tc::LineSearchTestCase) = minimum(tc.values)

function NLSolversBase.OnceDifferentiable(tc::LineSearchTestCase)
    function fdf(x)
        x < zero(x) && throw(ArgumentError("x must be nonnegative, got $x"))
        i = findfirst(>(x), tc.alphas)
        (i === nothing || x > tc.alphas[end-1]) && throw(ArgumentError("x must be <=$(tc.alphas[end-1]), got $x"))
        i -= 1
        xk = tc.alphas[i]
        xkp1 = tc.alphas[i + 1]
        dx = xkp1 - xk
        t = (x - xk) / dx
        h00t = 2t^3 - 3t^2 + 1
        h10t = t * (1 - t)^2
        h01t = t^2 * (3 - 2t)
        h11t = t^2 * (t - 1)
        val =
            h00t * tc.values[i] +
            h10t * dx * tc.slopes[i] +
            h01t * tc.values[i + 1] +
            h11t * dx * tc.slopes[i + 1]
        h00tp = 6t^2 - 6t
        h10tp = 3t^2 - 4t + 1
        h01tp = -6t^2 + 6 * t
        h11tp = 3t^2 - 2t
        slope =
            (
                h00tp * tc.values[i] +
                h10tp * dx * tc.slopes[i] +
                h01tp * tc.values[i + 1] +
                h11tp * dx * tc.slopes[i + 1]
            ) / dx
        return val, slope
    end
    return OnceDifferentiable(α -> fdf(α)[1], α -> fdf(α)[2], fdf, [0.0])
end

end