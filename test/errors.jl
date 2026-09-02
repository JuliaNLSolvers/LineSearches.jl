@testset "Invalid parameters" begin
    @test_throws "ArgumentError: The Armijo constant must satisfy 0 < c_1 < 1. Got c_1 = 0.0" BackTracking(; c_1 = 0.0)
    @test_throws "ArgumentError: The Armijo constant must satisfy 0 < c_1 < 1. Got c_1 = 1.0" BackTracking(; c_1 = 1.0)
    @test_throws "ArgumentError: The backtracking factors must satisfy 0 < ρ_lo <= ρ_hi < 1. Got ρ_lo = 0.9 and ρ_hi = 0.5" BackTracking(; ρ_lo = 0.9, ρ_hi = 0.5)
    @test_throws "ArgumentError: The backtracking factors must satisfy 0 < ρ_lo <= ρ_hi < 1. Got ρ_lo = 0.1 and ρ_hi = 1.0" BackTracking(; ρ_hi = 1.0)
    @test_throws "ArgumentError: The interpolation order must be 2 or 3. Got order = 4" BackTracking(; order = 4)
    @test_throws "ArgumentError: The iteration limit must be positive. Got iterations = 0" BackTracking(; iterations = 0)
    @test_throws "ArgumentError: The maximum step length must be positive. Got maxstep = 0.0" BackTracking(; maxstep = 0.0)

    @test_throws "ArgumentError: The Wolfe constants must satisfy 0 < c_1 < c_2 < 1. Got c_1 = 0.9 and c_2 = 0.1" StrongWolfe(; c_1 = 0.9, c_2 = 0.1)
    @test_throws "ArgumentError: The Wolfe constants must satisfy 0 < c_1 < c_2 < 1. Got c_1 = 0.0001 and c_2 = 1.0" StrongWolfe(; c_2 = 1.0)
    # ρ = 1 leaves the bracket the same width forever
    @test_throws "ArgumentError: The bracket growth factor must be larger than one, otherwise bracketing never terminates. Got ρ = 1.0" StrongWolfe(; ρ = 1.0)
end

@testset "Invalid endpoints" begin
    ϕ(α) = (α - 1.0)^2
    dϕ(α) = 2 * (α - 1.0)
    ϕdϕ(α) = (ϕ(α), dϕ(α))

    @testset "$(nameof(typeof(ls)))" for ls in (BackTracking(), StrongWolfe(),
                                                MoreThuente(), HagerZhang())
        msg = "LineSearchException: Value and slope at step length = 0 must be finite. Step length: 0.0"
        @test_throws msg ls(ϕ, dϕ, ϕdϕ, 1.0, NaN, -2.0)
        @test_throws msg ls(ϕ, dϕ, ϕdϕ, 1.0, 1.0, NaN)

        # An ascent direction can never satisfy sufficient decrease
        @test_throws "LineSearchException: Search direction is not a direction of descent. Step length: 0.0" ls(
            α -> 1.0 + α, α -> 1.0, α -> (1.0 + α, 1.0), 1.0, 1.0, 1.0)
    end

    # BackTracking still accepts a zero slope, where Armijo is a plain decrease test
    @test BackTracking()(ϕ, dϕ, ϕdϕ, 1.0, ϕ(0.0), 0.0) == (1.0, ϕ(1.0))
end

@testset "cstep outside its bracket" begin
    # An internal invariant, so it reports like the rest of the package
    @test_throws "LineSearchException: MoreThuente: cstep called with alpha = 3.0 outside the bracket [0.0, 2.0]. Step length: 3.0" LineSearches.cstep(
        0.0, 1.0, -1.0, 2.0, 0.5, 0.2, 3.0, 0.4, 0.1, true, 0.0, 10.0)
end

# HagerZhang asserts alphas[iB] > alphas[iA] on what secant2! returns, so update! must
# not narrow a bracket onto one of its own endpoints. It picks between narrowing to c
# and keeping the bracket using a fresh evaluation at c, which can disagree with the
# stored slope or value there.
@testset "update! keeps a bracket of positive width" begin
    unused(α) = error("ϕdϕ must not be called on these branches")
    @testset "c = $(alphas[3])" for (alphas, values, slopes) in (
        ([0.0, 1.0, 0.0], [0.0, 0.5, 0.0], [-1.0, 0.5, 0.25]),      # on the lower end
        ([0.0, 1.0, 1.0], [0.0, 2.0, -0.5], [-1.0, -0.25, -0.25]),  # on the upper end
        ([0.0, 1.0, NaN], [0.0, 2.0, -0.5], [-1.0, -0.25, -0.25]),  # not finite
    )
        iA, iB, _ = LineSearches.update!(unused, alphas, values, slopes, 1, 2, 3, 0.0,
                                         0.0, -1.0, LineSearches.DEFAULTDELTA,
                                         LineSearches.DEFAULTSIGMA)
        @test (alphas[iA], alphas[iB]) == (0.0, 1.0)
    end
end

# The step length these report is whatever the halving loop reached, so match a prefix
@testset "No finite value" begin
    @test_throws r"^LineSearchException: Static: failed to achieve finite new evaluation point\." Static()(α -> NaN, 1.0)
    @test_throws r"^LineSearchException: Backtracking: failed to achieve finite new evaluation point\." BackTracking()(α -> NaN, 1.0, 1.0, -1.0)
end
