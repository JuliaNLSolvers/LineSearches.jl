mutable struct LineSearchException{T<:Real} <: Exception
    message::AbstractString
    alpha::T
end

abstract type AbstractLineSearch end

# For debugging
struct LineSearchCache{T}
    alphas::Vector{T}
    values::Vector{T}
    slopes::Vector{T}
end
"""
    cache = LineSearchCache{T}()

Initialize an empty cache for storing intermediate results during line search.
The `α`, `ϕ(α)`, and possibly `dϕ(α)` values computed during line search are
available in `cache.alphas`, `cache.values`, and `cache.slopes`, respectively.

# Example

```jldoctest
julia> ϕ(x) = (x - π)^4; dϕ(x) = 4*(x-π)^3;

julia> cache = LineSearchCache{Float64}();

julia> ls = BackTracking(; cache);

julia> ls(ϕ, 10.0, ϕ(0), dϕ(0))
(1.8481462933284658, 2.7989406670901373)

julia> cache
LineSearchCache{Float64}([0.0, 10.0, 1.8481462933284658], [97.40909103400242, 2212.550050116452, 2.7989406670901373], [-124.02510672119926])
```

Because `BackTracking` doesn't use derivatives except at `α=0`, only the initial slope was stored in the cache.
Other methods may store all three.
"""
LineSearchCache{T}() where T = LineSearchCache{T}(T[], T[], T[])
