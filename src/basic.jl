# This linesearch does not perform a search at all, but takes the step
# alpha in the search direction. By default, alpha = 1.0 is the full step.
# This algorithm is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.

function basic!{T}(df,
                   x::Vector{T},
                   s::Vector,
                   x_scratch::Vector,
                   gr_scratch::Vector,
                   lsr::LineSearchResults,
                   alpha::Real = 1.0,
                   mayterminate::Bool = false)
    @assert alpha > 0
    push!(lsr.alpha, alpha)

    # Count number of parameters
    n = length(x)
    # Move a distance of alpha in the direction of s
    @simd for i in 1:n
        @inbounds x_scratch[i] = x[i] + alpha * s[i]
    end

    # Evaluate f(x) at new position
    f_x_scratch = df.f(x_scratch)
    push!(lsr.value, f_x_scratch)
    f_calls = 1

    g_calls = 0

    return alpha, f_calls, g_calls
end
