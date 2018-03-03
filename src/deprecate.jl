# Deprecation warnings

# alphatry and alphainit has been replaced with InitialHagerZhang

const dep_alphatry = Ref(false)
const dep_alphainit = Ref(false)

function alphainit(alpha::Real,
                   x::AbstractArray{T},
                   gr::AbstractArray,
                   f_x::Real,
                   psi0::T = convert(T,0.01)) where T
    if dep_alphainit[] == false
      warn("`alphainit` is deprecated, use `InitialHagerZhang()` instead")
      dep_alphainit[] = true
    end
    if isnan(alpha)
        alpha = LineSearches._hzI0(x,gr,f_x,psi0)
    end
    alpha
end

function alphatry(alpha::T,
                  df,
                  x::AbstractArray,
                  s::AbstractArray,
                  xtmp::AbstractArray,
                  lsr::LineSearchResults,
                  psi1::Real = convert(T,0.2),
                  psi2::Real = convert(T,2),
                  psi3::Real = convert(T,0.1),
                  iterfinitemax::Integer = ceil(Integer, -log2(eps(T))),
                  alphamax::Real = convert(T, Inf),
                  verbose::Bool = false) where T
   if dep_alphatry[] == false
      warn("`alphatry` is deprecated, use `InitialHagerZhang()` instead")
      dep_alphatry[] = true
    end

    LineSearches._hzI12(alpha, df, x, s, xtmp, lsr, psi1, psi2, psi3, alphamax, verbose)
end


(ls::Static)(df, x, s, x_scratch, lsr, alpha, mayterminate) = (ls::Static)(df, x, s, x_scratch, alpha)

_static!(df, x::AbstractArray{T}, s::AbstractArray{T}, lsr::LineSearchResults,
         x_scratch::AbstractArray{T}, alpha::Real = 1.0, mayterminate::Bool = false) where T =
_static!(df, x, s, x_scratch, alpha)

(ls::MoreThuente)(df, x, s, x_new, lsr::LineSearchResults, stp, mayterminate) = _morethuente!(df, x, s, x_new, lsr.value[1], lsr.slope[1], stp, mayterminate; f_tol=ls.f_tol, gtol=ls.gtol, x_tol=ls.x_tol, stpmin=ls.alphamin, stpmax=ls.alphamax, maxfev=ls.maxfev)

(ls::BackTracking)(df, x, s, x_scratch, lsr::LineSearchResults, alpha, mayterminate) = _backtracking!(df, x, s, x_scratch, lsr.value[1], lsr.slope[1], alpha, mayterminate, ls.c1, ls.rhohi, ls.rholo, ls.iterations, ls.order, ls.maxstep)


_strongwolfe!(df, x, p, x_new, lsr::LineSearchResults, alpha0, mayterminate; c1 = 1e-4, c2 = 0.9, rho = 2.0) = _strongwolfe!(df, x, p, x_new, lsr.value[1], lsr.slope[1], alpha0, mayterminate; c1 = 1e-4, c2 = 0.9, rho = 2.0)

_hagerzhang!(df, x, s, xtmp, lsr::LineSearchResults{T}, c, mayterminate, delta = DEFAULTDELTA,
                       sigma = DEFAULTSIGMA, alphamax = convert(T,Inf), rho = convert(T,5),
                       epsilon = convert(T,1e-6), gamma = convert(T,0.66),
                       linesearchmax = 50, psi3 = convert(T,0.1), display = 0) where T = _hagerzhang!(df, x, s, xtmp, lsr.value[1], lsr.slope[1],  c, mayterminate, delta, sigma, alphamax, rho ,epsilon,gamma ,linesearchmax,psi3,display)

# This one will have state.lsr
function (is::InitialHagerZhang)(state, dphi_0, df)
   if isnan(state.f_x_previous) && isnan(is.α0)
       # If we're at the first iteration (f_x_previous is NaN)
       # and the user has not provided an initial step size (is.α0 is NaN),
       # then we
       # pick the initial step size according to HZ #I0
       state.alpha = _hzI0(state.x, NLSolversBase.gradient(df),
                           NLSolversBase.value(df),
                           convert(eltype(state.x), is.ψ0)) # Hack to deal with type instability between is{T} and state.x
       state.mayterminate = false
   else
       # Pick the initial step size according to HZ #I1-2
       state.alpha, state.mayterminate =
           _hzI12(state.alpha, df, state.x, state.s, state.x_ls, state.lsr,
                  is.ψ1, is.ψ2, is.ψ3, is.αmax, is.verbose)
   end
   return state.alpha
end


_hzI12(alpha::T, df, x, s, xtmp, lsr::LineSearchResults{T}, psi1 = convert(T,0.2),
              psi2 = convert(T,2.0), psi3 = convert(T,0.1), alphamax = convert(T, Inf),
              verbose = false) where T = _hzI12(alpha, df, x, s, xtmp,
              lsr.value[1], lsr.slope[1], psi1, psi2, psi3, alphamax, verbose)
