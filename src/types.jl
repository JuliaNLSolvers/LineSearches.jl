# A cache for results from line search methods (to avoid recomputation)
type LineSearchResults{T}
    alpha::Vector{T}
    value::Vector{T}
    slope::Vector{T}
    nfailures::Int
end

LineSearchResults{T}(::Type{T}) = LineSearchResults(T[], T[], T[], 0)

Base.length(lsr::LineSearchResults) = length(lsr.alpha)

function Base.push!{T}(lsr::LineSearchResults{T}, a, v, d)
    push!(lsr.alpha, convert(T, a))
    push!(lsr.value, convert(T, v))
    push!(lsr.slope, convert(T, d))
    return
end

function clear!(lsr::LineSearchResults)
    empty!(lsr.alpha)
    empty!(lsr.value)
    empty!(lsr.slope)
    return
    # nfailures is deliberately not set to 0
end

type LineSearchException{T<:Real} <: Exception
    message::AbstractString
    alpha::T
    lsr::LineSearchResults
end
