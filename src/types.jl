mutable struct LineSearchException{T<:Real} <: Exception
    message::AbstractString
    alpha::T
end
