# A cache for results from line search methods (to avoid recomputation)
mutable struct LineSearchResults{T}
    alpha::Vector{T}
    value::Vector{T}
    slope::Vector{T}
    nfailures::Int
end

LineSearchResults(::Type{T}) where {T} = LineSearchResults(T[], T[], T[], 0)

Base.length(lsr::LineSearchResults) = length(lsr.alpha)

function Base.push!(lsr::LineSearchResults{T}, a, v, d) where T
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

mutable struct LineSearchException{T<:Real} <: Exception
    message::AbstractString
    alpha::T
    lsr::LineSearchResults
end
