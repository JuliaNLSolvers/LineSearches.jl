"""
`BackTracking` specifies a backtracking line-search that uses
a quadratic or cubic interpolant to determine the reduction in step-size.
E.g.,
if f(α) > f(0) + c₁ α f'(0), then the quadratic interpolant of
f(0), f'(0), f(α) has a minimiser α' in the open interval (0, α). More strongly,
there exists a factor ρ = ρ(c₁) such that α' ≦ ρ α.

This is a modification of the algorithm described in Nocedal Wright (2nd ed), Sec. 3.5.
"""
struct BackTracking{TF,TI}
  c_1::TF
  ρ_hi::TF
  ρ_lo::TF
  iterations::TI
  order::TI
  maxstep::TF
end

function BackTracking(args...; kwargs...)
  BackTracking{Float64}(args...; kwargs...)
end

function BackTracking{TF}(args...; kwargs...) where {TF}
  BackTracking{TF,Int}(args...; kwargs...)
end

function BackTracking{TF,TI}(;
  c_1::Real=1.0e-4, ρ_hi::Real=0.5, ρ_lo::Real=0.1,
  iterations::Integer=1_000, order::Integer=3, maxstep::Real=Inf
) where {TF,TI}
  # Impose 0 < ρ_hi < 1
  if (ρ_hi <= 0) || (ρ_hi >= 1)
    ρ_hi = 1 / 2
    msg = """
    The upper bound for the backtracking factor has to lie between
    0 and 1.
    Setting ρ_hi = $(ρ_hi).
    """
    warn(msg)
  end

  # Impose 0 < ρ_lo <= ρ_hi < 1
  if (ρ_lo <= 0) || (ρ_lo >= 1) || (ρ_lo > ρ_hi)
    ρ_lo = ρ_hi / 5
    msg = """
    The lower bound for the backtracking factor has to lie between
    0 and 1, and be smaller than the upper bound.
    Setting ρ_lo = $(ρ_lo).
    """
    warn(msg)
  end

  # Impose positive number of maximum iterations
  if (iterations <= 0) || !isinteger(iterations)
    iterations = trunc(Int, iterations)
    if iterations <= 0
      iterations = 1_000
    end
    msg = """
    The number of maximum iterations has to be a positive integer.
    Setting iterations = $(iterations).
    """
    warn(msg)
  end

  # Impose order in (2, 3)
  if order != 2 && order != 3
    order = 3
    msg = """The order has to be either 2 or 3.
    Setting order = $(order).
    """
    warn(msg)
  end

  # Impose maxstep
  if !isreal(maxstep)
    maxstep = Inf
    msg = """The maximum step size has to be real.
    Setting maxstep = $(maxstep).
    """
    warn(msg)
  end

  # Impose c_1 > 0
  if c_1 < 0
    c_1 = 1.0e-4
    msg = "The Armijo constant hast to be positive.
    Setting c_1 = $(c_1)."
    warn(msg)
  end

  # # Impose backtracking factor (for order = 2)
  # # The quadratic update rule come with a backtracking factor
  # #    ρ = 1 / 2 / (1 - c_1).
  # # We need c_1 > 0, and we want ρ < 1,
  # # so 0 < c_1 < 1/2 and 1/2 < ρ < 1.
  # ρ = ρ_hi # Could take another choice here
  # c_1_ρ = 1 - 1 / (2 * ρ)
  # if c_1 > c_1_ρ
  #   c_1 = c_1_ρ
  #   msg = """The Armijo constant c_1 is too large.
  #   Setting c_1 = $(c_1_ρ)."""
  #   warn(msg)
  # end

  BackTracking{TF,TI}(c_1, ρ_hi, ρ_lo, iterations, order, maxstep)
end

function (ls::BackTracking)(
  df::AbstractObjective, x::AbstractArray{T}, s::AbstractArray{T},
  α_0::Tα=real(T)(1), x_new::AbstractArray{T}=similar(x),
  ϕ_0=nothing, dϕ_0=nothing, alphamax=convert(real(T), Inf)
) where {T,Tα}
  ϕ, dϕ = make_ϕ_dϕ(df, x_new, x, s)

  if isnothing(ϕ_0)
    ϕ_0 = ϕ(Tα(0))
  end
  if isnothing(dϕ_0)
    dϕ_0 = dϕ(Tα(0))
  end

  α_0 = min(α_0, min(alphamax, ls.maxstep / norm(s, Inf)))
  ls(ϕ, α_0, ϕ_0, dϕ_0)
end

function (ls::BackTracking)(ϕ, dϕ, ϕdϕ, α_0, ϕ_0, dϕ_0)
  ls(ϕ, α_0, ϕ_0, dϕ_0)
end

# TODO: Should we deprecate the interface that only uses the ϕ argument?
function (ls::BackTracking)(ϕ, α_0::Tα, ϕ_0, dϕ_0) where {Tα}
  @unpack c_1, ρ_hi, ρ_lo, iterations, order = ls
  ε = eps(real(Tα))

  # Initialise α_1 and α_2
  α_1, α_2 = α_0, α_0
  ϕ_1, ϕ_2 = ϕ_0, ϕ_0

  # Backtrack until ϕ(α_2) is finite
  iterfinite = 0
  iterfinitemax = -log2(ε)
  ϕ_2 = ϕ(α_1)
  while !isfinite(ϕ_2) && iterfinite < iterfinitemax
    iterfinite += 1
    α_1 = α_2
    α_2 = α_1 / 2
    ϕ_2 = ϕ(α_2)
  end

  # Backtrack until sufficient decrease
  iteration = 0
  while (ϕ_2 > ϕ_0 + c_1 * α_2 * dϕ_0) && (iteration <= iterations)
    iteration += 1

    # Shrink proposed step-size:
    if (order == 2) || (iteration == 1)
      # Backtracking via quadratic interpolation:
      # This interpolates the available data
      #    ϕ(0), ϕ'(0), ϕ(α)
      # with a quadractic which is then minimised.
      α_tmp = -(dϕ_0 * α_2^2) / (2 * (ϕ_2 - ϕ_0 - dϕ_0 * α_2))
    else
      # Backtracking via cubic interpolation:
      # This interpolates the available data
      #    ϕ(0), ϕ'(0), ϕ(α_1), ϕ(α_2)
      # with a cubic function which is then minimised.
      α_1², α_2² = α_1^2, α_2^2
      α_1³, α_2³ = α_1² * α_1, α_2² * α_2

      δ_1 = ϕ_1 - ϕ_0 - dϕ_0 * α_1
      δ_2 = ϕ_2 - ϕ_0 - dϕ_0 * α_2

      invdet = one(Tα) / (α_1² * α_2² * (α_2 - α_1))

      a = (α_1² * δ_2 - α_2² * δ_1) * invdet
      b = (α_2³ * δ_1 - α_1³ * δ_2) * invdet

      if isapprox(a, zero(Tα), atol=ε)
        # Degenerate quadratic case
        α_tmp = -dϕ_0 / (2 * b)
      else
        # General cubic case, avoiding numerical cancellation
        Δ = max(b^2 - 3 * a * dϕ_0, zero(Tα))
        α_tmp = -dϕ_0 / (b + sqrt(Δ))
      end
    end

    # Clamp α_tmp to avoid too small / large reductions
    α_tmp = NaNMath.min(α_tmp, α_2 * ρ_hi)
    α_tmp = NaNMath.max(α_tmp, α_2 * ρ_lo)

    # Update (α_1, α_2)
    α_1, α_2 = α_2, α_tmp
    ϕ_1, ϕ_2 = ϕ_2, ϕ(α_2)
  end

  # Ensure termination
  if iteration > iterations
    msg = "Linesearch failed to converge, reached maximum iterations $(iterations)."
    throw(LineSearchException(msg, α_2))
  end

  return α_2, ϕ_2
end
