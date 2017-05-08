# Deprecation warnings

bt3!{T}(df,
        x::Vector{T},
        s::Vector,
        x_scratch::Vector,
        gr_scratch::Vector,
        lsr::LineSearchResults,
        alpha::Real = 1.0,
        mayterminate::Bool = false,
        c1::Real = 1e-4,
        rhohi::Real = 0.5,
        rholo::Real = 0.1,
        iterations::Integer = 1_000) =
            backtracking!(df,x,s,x_scratch,gr_scratch,
                          lsr,alpha,mayterminate,c1,
                          rhohi,rholo,iterations,3)

bt2!{T}(df,
        x::Vector{T},
        s::Vector,
        x_scratch::Vector,
        gr_scratch::Vector,
        lsr::LineSearchResults,
        alpha::Real = 1.0,
        mayterminate::Bool = false,
        c1::Real = 1e-4,
        rhohi::Real = 0.5,
        rholo::Real = 0.1,
        iterations::Integer = 1_000) =
            backtracking!(df,x,s,x_scratch,gr_scratch,
                          lsr,alpha,mayterminate,c1,
                          rhohi,rholo,iterations,2)
