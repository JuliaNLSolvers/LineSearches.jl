# Deprecation warnings

# >>>>>>>> Deprecation warnings for linesearch functions (vs types)

export hagerzhang!, backtracking!, strongwolfe!,
    morethuente!, bt2!, bt3!, basic!

const dep_bt2 = Ref(false)
const dep_bt3 = Ref(false)
const dep_backtracking = Ref(false)
const dep_basic = Ref(false)
const dep_strongwolfe = Ref(false)
const dep_morethuente = Ref(false)
const dep_hagerzhang = Ref(false)

function _deprecate(msg, depvar, lsfun, args...)
   if depvar[] == false
      warn(msg)
      depvar[] = true
   end
   return lsfun(args...)
end


_bt3!(args...) = _backtracking!(args...)

_bt2!(df, x, s, x_scratch, lsr, alpha, mayterminate::Bool,
      c1::Real = 1e-4, rhohi::Real = 0.5, rholo::Real = 0.1, iterations::Integer = 1_000) =
          _backtracking!(df, x, s, x_scratch, lsr, alpha, mayterminate,
                         c1, rhohi, rholo, iterations, 2)


bt3!(args...) = _deprecate(
   "`bt3!` is deprecated, use `BackTracking()` instead",
   dep_bt3, _bt3!, args...)

bt2!(args...) = _deprecate(
   "`bt2!` is deprecated, use `BackTracking(order = 2)` instead",
   dep_bt2, _bt2!, args...)

backtracking!(args...) = _deprecate(
   "`backtracking!` is deprecated, use `BackTracking()` instead",
   dep_backtracking, _backtracking!, args...)

strongwolfe!(args...) = _deprecate(
   "`strongwolfe!` is deprecated, use `StrongWolfe()` instead",
   dep_strongwolfe, _strongwolfe!, args...)

morethuente!(args...) = _deprecate(
   "`morethuente!` is deprecated, use `MoreThuente()` instead",
   dep_morethuente, _morethuente!, args...)

hagerzhang!(args...) = _deprecate(
   "`hagerzhang!` is deprecated, use `HagerZhang()` instead",
   dep_hagerzhang, _hagerzhang!, args...)

basic!(args...) = _deprecate(
   "`basic!` is deprecated, use `Static()` instead",
   dep_basic, _static!, args...)

# <<<<<<<<<<<< end deprecation warnings for linesearch functions


# <<<< Start deprecation of gradient storage removal

const dep_g_bt2 = Ref(false)
const dep_g_bt3 = Ref(false)
const dep_g_backtracking = Ref(false)
const dep_g_static = Ref(false)
const dep_g_strongwolfe = Ref(false)
const dep_g_morethuente = Ref(false)
const dep_g_hagerzhang = Ref(false)
const dep_g_alphatry = Ref(false)

function _warn_g(depvar)
    if depvar[] == false
        warn("You no longer have to provide a 'g'(gradient storage) input")
        depvar[] = true
    end
end

function _hagerzhang!(df, x, s, xtmp, g, lsr, c, mayterminate, args...)
    _warn_g(dep_g_hagerzhang)
    retval = _hagerzhang!(df, x, s, xtmp, lsr, c, mayterminate, args...)
    copy!(g, df.g)
    return retval
end

function _backtracking!(df, x, s, xtmp, g, lsr, c, mayterminate, args...)
    _warn_g(dep_g_backtracking)
    retval = _backtracking!(df, x, s, xtmp, lsr, c, mayterminate, args...)
    copy!(g, df.g)
    return retval
end

function _bt2!(df, x, s, xtmp, g, lsr, c, mayterminate, args...)
    _warn_g(dep_g_bt2)
    retval = _bt2!(df, x, s, xtmp, lsr, c, mayterminate, args...)
    copy!(g, df.g)
    return retval
end

_bt2!(df, x, s, x_scratch, g, lsr, alpha, mayterminate,
      c1::Real = 1e-4, rhohi::Real = 0.5, rholo::Real = 0.1, iterations::Integer = 1_000) =
          _backtracking!(df, x, s, x_scratch, g, lsr, alpha, mayterminate,
                         c1, rhohi, rholo, iterations, 2)

function _bt3!(df, x, s, xtmp, g, lsr, c, mayterminate, args...)
    _warn_g(dep_g_bt3)
    retval = _bt3!(df, x, s, xtmp, lsr, c, mayterminate, args...)
    copy!(g, df.g)
    return retval
end

function _static!(df, x, s, xtmp, g, lsr, c, mayterminate, args...)
    _warn_g(dep_g_static)
    retval = _static!(df, x, s, xtmp, lsr, c, mayterminate, args...)
    copy!(g, df.g)
    return retval
end

function _strongwolfe!(df, x, s, xtmp, g, lsr, c, mayterminate, args...)
    _warn_g(dep_g_strongwolfe)
    retval = _strongwolfe!(df, x, s, xtmp, lsr, c, mayterminate, args...)
    copy!(g, df.g)
    return retval
end

function _morethuente!(df, x, s, xtmp, g, lsr, c, mayterminate, args...)
    _warn_g(dep_g_morethuente)
    retval = _morethuente!(df, x, s, xtmp, lsr, c, mayterminate, args...)
    copy!(g, df.g)
    return retval
end

function alphatry(alpha::T, df, x::Array, s::Array, xtmp::Array, g::Array, lsr::LineSearchResults, args...) where T
    _warn_g(dep_g_alphatry)
    alphatry(alpha, df, x, s, xtmp, lsr, args...)
end

# >>>> End deprecation  f gradient storage removal
