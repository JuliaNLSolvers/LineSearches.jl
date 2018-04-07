@with_kw struct LineSearch{LS, IS}
    linesearch!::LS = BackTracking()
    initial_step::IS = InitialStatic()
end
function (ls::LineSearch)(ϕ, dϕ, α_0 = 1.0,
                          ϕ_0 = ϕ(α_0), dϕ_0 = dϕ(α_0))
    @unpack linesearch!, initial_step = ls
    if display
        println("New linesearch")
    end
    alpha > T(0) || error("Initial step-size must be positive")
    (isfinite(phi_0) && isfinite(dphi_0)) || error("Initial value and slope must be finite")

end
