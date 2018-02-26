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


(ls::MoreThuente)(args...) =
       _morethuente!(args...;
                   f_tol=ls.f_tol, gtol=ls.gtol, x_tol=ls.x_tol, stpmin=ls.alphamin,
                   stpmax=ls.alphamax, maxfev=ls.maxfev)
(ls::MoreThuente)(df, x, s, x_new, lsr::LineSearchResults, stp, mayterminate) =
       _morethuente!(df, x, s, x_new, lsr.value[1], lsr.slope[1], stp, mayterminate;
                   f_tol=ls.f_tol, gtol=ls.gtol, x_tol=ls.x_tol, stpmin=ls.alphamin,
                   stpmax=ls.alphamax, maxfev=ls.maxfev)


(ls::BackTracking)(df, x, s, x_scratch, lsr::LineSearchResults, alpha, mayterminate) =
   _backtracking!(df, x, s, x_scratch, lsr.value[1], lsr.slope[1], alpha, mayterminate,
            ls.c1, ls.rhohi, ls.rholo, ls.iterations, ls.order, ls.maxstep)


_strongwolfe!(df, x::AbstractArray{T}, p::AbstractArray{T}, x_new::AbstractArray{T},
                      lsr::LineSearchResults{T}, alpha0::Real, mayterminate::Bool;
                      c1::Real = 1e-4, c2::Real = 0.9, rho::Real = 2.0) where T =
                      _strongwolfe!(df, x, p, x_new, lsr.value[1], lsr.slope[1], alpha0, mayterminate; c1 = 1e-4, c2 = 0.9, rho = 2.0)

_hagerzhang!(df, x::AbstractArray{T}, s::AbstractArray{T}, xtmp::AbstractArray{T}, lsr::LineSearchResults, c::Real, mayterminate::Bool, delta::Real = DEFAULTDELTA,
                       sigma::Real = DEFAULTSIGMA, alphamax::Real = convert(T,Inf), rho::Real = convert(T,5),
                       epsilon::Real = convert(T,1e-6), gamma::Real = convert(T,0.66),
                       linesearchmax::Integer = 50, psi3::Real = convert(T,0.1), display::Integer = 0) where T =
                       _hagerzhang!(df, x, s, xtmp, lsr.value[1], lsr.slope[1],  c, mayterminate, delta, sigma, alphamax, rho ,epsilon,gamma ,linesearchmax,psi3,display)

_hzI12(alpha::T, df, x::AbstractArray{T}, s::AbstractArray{T},
              xtmp::AbstractArray{T},
              lsr::LineSearchResults,
              psi1::Real = convert(T,0.2),
              psi2::Real = convert(T,2.0),
              psi3::Real = convert(T,0.1),
              alphamax::Real = convert(T, Inf),
              verbose::Bool = false) where T = _hzI12(alpha,
              df,
              x,
              s,
              xtmp,
              lsr.value[1],
              lsr.slope[1],
              psi1,
              psi2,
              psi3,
              alphamax,
              verbose)
