@with_kw struct LineSearch{LS, IS}
    linesearch!::LS = BackTracking()
    initial_step::IS = InitialStatic()
end
function (ls::LineSearch)(ϕ, dϕ,
                          α_0 = 1.0,
                          ϕ_0 = ϕ(α_0), dϕ_0 = dϕ(α_0);
                          step = nothing, maxstep = nothing,
                          maxalpha = Inf, display = false)
    # Let the user know that we've entered the line search
    if display
      println("* New line search")
    end

    # Error handling

    # Make sure maxstep was input correctly (neither can be omitted if the other
    # one present)
@compat    if !(typeof(step) <: Void) || !(typeof(maxstep) <: Void)
        if step <: Void
            throw(ErrorException("When you set the maxstep keyword, you also need to provide the actual step using the `step` keyword."))
        elseif maxstep <: Void
            throw(ErrorException("When you set the step keyword, you also need to provide the largest allowed step using the `maxstep` keyword."))
        end
        maxalpha = norm(step, Inf)/maxstep
    end

    # Make sure that the initial step size is positive
    α_0 > T(0) || error("Initial step size must be positive")

    # Make sure that initial step size, value and slope are finite
    all(isfinite.(α_0, phi_0, dphi_0)) || error("Initial step size, value and slope must be finite")

    # Error handling over; unpack ls
    @unpack linesearch!, initial_step = ls


end
