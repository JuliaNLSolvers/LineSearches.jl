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
