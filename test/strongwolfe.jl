@testset "interpolate" begin
    interp = LineSearches.interpolate

    # Each of these used to throw, return NaN, or land outside (0.5, 2.0)
    @test interp(0.5, 2.0, 1.0, 0.0, -1.0, -1.0) == 1.25      # negative radicand
    @test interp(0.5, 2.0, 1.0, 0.2, -1e160, 1e160) == 1.25   # d1^2 overflows
    @test interp(0.5, 2.0, 1.0, 0.2, -1e-180, 1e-180) == 1.25 # slope product underflows
    @test interp(0.5, 2.0, 0.0, 0.0, 0.0, 0.0) == 1.25        # nothing to interpolate

    # Where the cubic is well posed the safeguards cost no accuracy
    args = (0.0, 1.0, 1.0, 0.2, -3.0, 0.7)
    @test interp(args...) ≈ Float64(interp(big.(args)...)) rtol = 1e-15
end

@testset "zoom with same-sign slopes" begin
    # The first zoom call site brackets on ϕ alone, so both slopes may share a sign
    ϕ(x) = (x - π)^4
    dϕ(x) = 4 * (x - π)^3

    a, val = LineSearches.zoom(0.0, ϕ(0.0), dϕ(0.0), 1.0, ϕ(1.0), dϕ(1.0),
                               dϕ(0.0), ϕ(0.0), x -> (ϕ(x), dϕ(x)), nothing, 1e-4, 0.9)
    @test 0 < a < 1
    @test val == ϕ(a)
    @test abs(dϕ(a)) <= 0.9 * abs(dϕ(0.0))
end
