# TODO: Implement safeguards


"""
`StrongWolfe`: This linesearch algorithm guarantees that the step length
satisfies the (strong) Wolfe conditions.
See Nocedal and Wright - Algorithms 3.5 and 3.6

This algorithm is mostly of theoretical interest, users should most likely
use `MoreThuente`, `HagerZhang` or `BackTracking`.

## Parameters:  (and defaults)
* `c_1 = 1e-4`: Armijo condition
* `c_2 = 0.9` : second (strong) Wolfe condition
* `ρ = 2.0` : bracket growth
"""
@with_kw struct StrongWolfe{T} <: AbstractLineSearch
	c_1::T = 1e-4
	c_2::T = 0.9
	ρ::T = 2.0
	cache::Union{Nothing, LineSearchCache{T}} = nothing
end

"""
	(ls::StrongWolfe)(df, x::AbstractArray, p::AbstractArray, alpha0::Real, x_new, ϕ_0, dϕ_0) -> alpha, ϕalpha

Given a differentiable function `df` (in the sense of `NLSolversBase.OnceDifferentiable` or
`NLSolversBase.TwiceDifferentiable`), a multidimensional starting point `x` and step `p`,
and a guess `alpha0` for the step length, find an `alpha` satisfying the strong Wolfe conditions.

See the one-dimensional method for additional details.
"""
function (ls::StrongWolfe)(df, x::AbstractArray{T},
	p::AbstractArray{T}, α::Real, x_new::AbstractArray{T},
	ϕ_0, dϕ_0) where T
	ϕ, dϕ, ϕdϕ = make_ϕ_dϕ_ϕdϕ(df, x_new, x, p)
	ls(ϕ, dϕ, ϕdϕ, α, ϕ_0, dϕ_0)
end

"""
	(ls::StrongWolfe)(ϕ, dϕ, ϕdϕ, alpha0, ϕ_0, dϕ_0) -> alpha, ϕalpha

Given `ϕ(alpha::Real)`, its derivative `dϕ`, a combined-evaluation function
`ϕdϕ(alpha) -> (ϕ(alpha), dϕ(alpha))`, and an initial guess `alpha0`,
identify a value of `alpha > 0` satisfying the strong Wolfe conditions.

`ϕ_0` and `dϕ_0` are the value and derivative, respectively, of `ϕ` at `alpha = 0.`

Both `alpha` and `ϕ(alpha)` are returned.
"""
function (ls::StrongWolfe)(ϕ, dϕ, ϕdϕ,
	alpha0::T, ϕ_0, dϕ_0) where T <: Real
	@unpack c_1, c_2, ρ, cache = ls
	emptycache!(cache)

	pushcache!(cache, zero(T), ϕ_0, dϕ_0)

	# Step-sizes
	a_0 = zero(T)
	a_iminus1 = a_0
	a_i = alpha0
	a_max = convert(T, 1000.0)

	# ϕ(alpha) = df.f(x + alpha * p)
	ϕ_a_iminus1 = ϕ_0
	ϕ_a_i = convert(T, NaN)

	# ϕ'(alpha) = dot(g(x + alpha * p), p)
	dϕ_a_i = convert(T, NaN)
	i = 1
	while a_i < a_max
		ϕ_a_i = ϕ(a_i)

		# Test Wolfe conditions
		if (ϕ_a_i > ϕ_0 + c_1 * a_i * dϕ_0) ||
		   (ϕ_a_i >= ϕ_a_iminus1 && i > 1)
			a_star = zoom(a_iminus1, a_i,
				dϕ_0, ϕ_0,
				ϕ, dϕ, ϕdϕ, cache)
			return a_star, ϕ(a_star)
		end

		dϕ_a_i = dϕ(a_i)
		pushcache!(cache, a_i, ϕ_a_i, dϕ_a_i)

		# Check condition 2
		if abs(dϕ_a_i) <= -c_2 * dϕ_0
			return a_i, ϕ_a_i
		end

		# Check condition 3
		if dϕ_a_i >= zero(T) # FIXME untested!
			a_star = zoom(a_i, a_iminus1,
				dϕ_0, ϕ_0, ϕ, dϕ, ϕdϕ, cache)
			return a_star, ϕ(a_star)
		end

		# Choose a_iplus1 from the interval (a_i, a_max)
		a_iminus1 = a_i
		a_i *= ρ

		# Update ϕ_a_iminus1
		ϕ_a_iminus1 = ϕ_a_i
		i += 1
	end

	error("StrongWolfe did not converge")
end

function zoom(
	a_lo::T,
	a_hi::T,
	dϕ_0::Real,
	ϕ_0::Real,
	ϕ,
	dϕ,
	ϕdϕ,
	cache,
	c_1::Real = 1.0 / 10.0^4,
	c_2::Real = 0.9) where T

	# Step-size
	a_j = convert(T, NaN)

	# Count iterations
	iteration = 0
	max_iterations = 10

	# Shrink bracket
	while iteration < max_iterations
		iteration += 1

		ϕ_a_lo, ϕprime_a_lo = ϕdϕ(a_lo)
		pushcache!(cache, a_lo, ϕ_a_lo, ϕprime_a_lo)

		ϕ_a_hi, ϕprime_a_hi = ϕdϕ(a_hi)

		pushcache!(cache, a_hi, ϕ_a_hi, ϕprime_a_hi)
		# Interpolate a_j
		if a_lo < a_hi
			a_j = interpolate(a_lo, a_hi,
				ϕ_a_lo, ϕ_a_hi,
				ϕprime_a_lo, ϕprime_a_hi)
		else
			# TODO: Check if this is needed
			a_j = interpolate(a_hi, a_lo,
				ϕ_a_hi, ϕ_a_lo,
				ϕprime_a_hi, ϕprime_a_lo)
		end

		# Evaluate ϕ(a_j)
		ϕ_a_j = ϕ(a_j)
		pushcache!(cache, a_j, ϕ_a_j)

		# Check Armijo
		if (ϕ_a_j > ϕ_0 + c_1 * a_j * dϕ_0) ||
		   (ϕ_a_j > ϕ_a_lo)
			a_hi = a_j
		else
			# Evaluate ϕprime(a_j)
			ϕprime_a_j = dϕ(a_j)
			if cache !== nothing
				push!(cache.slopes, ϕprime_a_j)
			end

			if abs(ϕprime_a_j) <= -c_2 * dϕ_0
				return a_j
			end

			if ϕprime_a_j * (a_hi - a_lo) >= zero(T)
				a_hi = a_lo
			end

			a_lo = a_j
		end
	end

	# Quasi-error response
	return a_j
end

# a_lo = a_{i - 1}
# a_hi = a_{i}
function interpolate(a_lo::Real, a_hi::Real,
	ϕ_a_lo::Real, ϕ_a_hi::Real,
	dϕ_a_lo::Real, dϕ_a_hi::Real)

	d1 = dϕ_a_lo + dϕ_a_hi - 3 * (ϕ_a_lo - ϕ_a_hi) / (a_lo - a_hi)
	d2 = sqrt(d1 * d1 - dϕ_a_lo * dϕ_a_hi)
	return a_hi - (a_hi - a_lo) *
				  ((dϕ_a_hi + d2 - d1) /
				   (dϕ_a_hi - dϕ_a_lo + 2 * d2))
end
