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

_bt2!(df, x, s, x_scratch, gr_scratch, lsr, alpha, mayterminate,
      c1::Real = 1e-4, rhohi::Real = 0.5, rholo::Real = 0.1, iterations::Integer = 1_000) =
      _backtracking!(df, x, s, x_scratch, gr_scratch, lsr, alpha, mayterminate,
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
