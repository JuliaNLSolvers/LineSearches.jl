# Standalone double-well comparison: HagerZhang vs HagerZhangLS
#
# φ(α) = (α² - 4)² - 5α
# φ'(α) = 4α(α² - 4) - 5
# φ(0) = 16, φ'(0) = -5 (descent direction ✓)

using Pkg
Pkg.activate(@__DIR__)
using LineSearches

# --- Test function ---
φdφ(α) = ((α^2 - 4)^2 - 5α, 4α * (α^2 - 4) - 5)
ϕ(α) = φdφ(α)[1]
dϕ(α) = φdφ(α)[2]
ϕdϕ(α) = φdφ(α)

φ0, dφ0 = φdφ(0.0)
c0 = 1.0

println("Double-well: φ(α) = (α²-4)² - 5α")
println("φ(0) = $φ0, φ'(0) = $dφ0, c₀ = $c0")
println()

# --- Run both with default σ = 0.9 ---
hz = HagerZhang()
hz.mayterminate[] = false
hzls = HagerZhangLS()

α_hz, fα_hz = hz(ϕ, dϕ, ϕdϕ, c0, φ0, dφ0)
α_hzls, fα_hzls = hzls(ϕ, dϕ, ϕdϕ, c0, φ0, dφ0)

println("=== Default (σ = 0.9) ===")
println("  HagerZhang:   α = $α_hz,  φ(α) = $fα_hz")
println("  HagerZhangLS: α = $α_hzls, φ(α) = $fα_hzls")
println()

# --- Run both with strict σ = 0.5 ---
hz_strict = HagerZhang(sigma = 0.5)
hz_strict.mayterminate[] = false
hzls_strict = HagerZhangLS(curvature = 0.5)

α_hz_s, fα_hz_s = hz_strict(ϕ, dϕ, ϕdϕ, c0, φ0, dφ0)
α_hzls_s, fα_hzls_s = hzls_strict(ϕ, dϕ, ϕdϕ, c0, φ0, dφ0)

println("=== Strict (σ = 0.5) ===")
println("  HagerZhang:   α = $α_hz_s,  φ(α) = $fα_hz_s")
println("  HagerZhangLS: α = $α_hzls_s, φ(α) = $fα_hzls_s")
println()

# --- Verify Wolfe conditions ---
function check_wolfe(α, φα, dφα, φ0, dφ0, δ, σ, ϵ)
    wolfe = α > 0 && δ * dφ0 ≥ (φα - φ0) / α && dφα ≥ σ * dφ0
    awolfe = (2δ - 1) * dφ0 ≥ dφα ≥ σ * dφ0 && φα ≤ φ0 + ϵ * abs(φ0)
    armijo = α > 0 && δ * dφ0 ≥ (φα - φ0) / α
    curvature = dφα ≥ σ * dφ0
    return (; wolfe, awolfe, armijo, curvature)
end

for (label, σ, αhz, αhzls) in [("σ=0.9", 0.9, α_hz, α_hzls), ("σ=0.5", 0.5, α_hz_s, α_hzls_s)]
    println("=== Wolfe check ($label) ===")
    for (name, α) in [("HZ", αhz), ("HZLS", αhzls)]
        if isfinite(α) && α > 0
            φα, dφα = φdφ(α)
            conds = check_wolfe(α, φα, dφα, φ0, dφ0, 0.1, σ, 1e-6)
            println("  $name (α=$α): wolfe=$(conds.wolfe), awolfe=$(conds.awolfe), armijo=$(conds.armijo), curvature=$(conds.curvature)")
        else
            println("  $name: failed (α=$α)")
        end
    end
    println()
end
