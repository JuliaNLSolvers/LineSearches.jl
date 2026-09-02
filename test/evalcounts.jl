@testset "Evaluation counts" begin
    # `only_fg!` passes `nothing` for the slot it does not want, so logging the requested
    # slots shows whether a fused or a split evaluation was used. NLSolversBase's own
    # f_calls/g_calls cannot tell them apart.
    f(x) = sum(abs2, x) + 0.1 * sum(z -> z^4, x)
    ∇f(x) = 2 .* x .+ 0.4 .* x .^ 3

    function counted()
        log = Symbol[]
        function fdf(F, G, x)
            push!(log, F === nothing ? :grad : (G === nothing ? :val : :fused))
            G === nothing || (G .= ∇f(x))
            F === nothing || return f(x)
            return nothing
        end
        return fdf, log
    end
    objective(fdf, x) =
        NLSolversBase.OnceDifferentiable(NLSolversBase.only_fg!(fdf), copy(x), 0.0)

    x = [1.0, 2.0, 3.0]
    s = -[2.0, 4.0, 6.0]
    ϕ_0 = f(x)
    dϕ_0 = sum(∇f(x) .* s)

    lss = (BackTracking(), BackTracking(; order = 2), LineSearches.HagerZhangLS())

    @testset "α=0: $(nameof(typeof(ls)))" for ls in lss
        fdf, log = counted()
        α, val = ls(objective(fdf, x), x, s, 1.0, similar(x))
        @test first(log) == :fused
        @test count(==(:grad), log) == 0
        @test α ≈ 0.5
        @test val ≈ f(x .+ α .* s)

        # A supplied endpoint keeps the narrower call
        fdf, log = counted()
        ls(objective(fdf, x), x, s, 1.0, similar(x), ϕ_0, nothing)
        @test first(log) == :grad

        fdf, log = counted()
        ls(objective(fdf, x), x, s, 1.0, similar(x), nothing, dϕ_0)
        @test first(log) == :val
    end

    # Steep descent with a flat tail, so a tight c_2 makes zoom iterate
    ψ(α) = -α / (1 + 10α^2) + 0.02α^2
    dψ(α) = -(1 - 10α^2) / (1 + 10α^2)^2 + 0.04α

    # zoom reassigns its endpoints only to already-evaluated points, so it costs 2
    # evaluations to seed the bracket plus 1 per iteration (was ~4 per iteration).
    @testset "zoom, c_2=$c_2" for (c_2, entries, expected) in
            ((0.5, 3, 0.412201859), (0.1, 4, 0.3461435202), (0.01, 7, 0.312411912))
        n = Ref(0)
        ϕdϕ = α -> (n[] += 1; (ψ(α), dψ(α)))
        a = LineSearches.zoom(0.0, 1.0, dψ(0.0), ψ(0.0), ϕdϕ, nothing, 1e-4, c_2)
        @test a ≈ expected
        @test n[] == entries
    end

    # zoom used to default c_1/c_2 in its signature and neither call site passed them,
    # so a tight c_2 was ignored and the returned step violated the curvature condition.
    @testset "c_2 is honoured, c_2=$c_2" for c_2 in (0.1, 0.05, 0.01)
        α, val = StrongWolfe(; c_2)(ψ, dψ, α -> (ψ(α), dψ(α)), 1.0, ψ(0.0), dψ(0.0))
        @test abs(dψ(α)) <= c_2 * abs(dψ(0.0))
        @test val ≈ ψ(α)
    end

    # HagerZhang stores every evaluated point, so it should never ask for the same α twice
    problems = ((ψ, dψ),
                (α -> (α^2 - 11)^2 + (α - 7)^2, α -> 4α * (α^2 - 11) + 2(α - 7)),
                (α -> 1 - 1 / (1 + (α - 2)^2), α -> 2(α - 2) / (1 + (α - 2)^2)^2))
    @testset "HagerZhang evaluates each α once" for (g, dg) in problems,
                                                    α0 in (0.1, 1.0, 5.0)
        seen = Float64[]
        α, val = HagerZhang()(α -> (push!(seen, α); (g(α), dg(α))), α0, g(0.0), dg(0.0))
        @test allunique(seen)
        @test val ≈ g(α)
    end
end
