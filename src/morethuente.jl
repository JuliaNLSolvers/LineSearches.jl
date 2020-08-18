#   Translation of Matlab version by John Myles White
#   Translation of minpack subroutine cvsrch
#   Dianne O'Leary   July 1991
#
#     **********
#
#     Subroutine cvsrch
#
#     The purpose of cvsrch is to find a step which satisfies
#     a sufficient decrease condition and a curvature condition.
#     The user must provide a subroutine which calculates the
#     function and the gradient.
#
#     At each stage the subroutine updates an interval of
#     uncertainty with endpoints stx and sty. The interval of
#     uncertainty is initially chosen so that it contains a
#     minimizer of the modified function
#
#          f(x + stp * s) - f(x) - f_tol * stp * (gradf(x)' * s).
#
#     If a step is obtained for which the modified function
#     has a nonpositive function value and nonnegative derivative,
#     then the interval of uncertainty is chosen so that it
#     contains a minimizer of f(x + stp * s).
#
#     The algorithm is designed to find a step which satisfies
#     the sufficient decrease condition
#
#           f(x + stp * s) <= f(x) + f_tol * stp * (gradf(x)' * s),
#
#     and the curvature condition
#
#           abs(gradf(x + stp * s)' * s)) <= gtol * abs(gradf(x)' * s).
#
#     If f_tol is less than gtol and if, for example, the function
#     is bounded below, then there is always a step which satisfies
#     both conditions. If no step can be found which satisfies both
#     conditions, then the algorithm usually stops when rounding
#     errors prevent further progress. In this case stp only
#     satisfies the sufficient decrease condition.
#
#     The subroutine statement is
#
#        subroutine cvsrch(df,n,x,f,s,stp,f_tol,gtol,x_tol,
#                          alphamin,alphamax,maxfev,info,nfev,wa)
#
#     where
#
# df is the name of the user-supplied subroutine which
#    calculates the function and the gradient.  df must
#    be declared in an external statement in the user
#    calling program, and should be written as follows.
#
#    function [f,g] = df(n,x) (Matlab)
#                     (10/2010 change in documentation)
#                     (derived from Fortran subroutine df(n,x,f,g))
#    integer n
#    f
#    x(n),g(n)
#
#    Calculate the function at x and
#    return this value in the variable f.
#    Calculate the gradient at x and
#    return this vector in g.
#
#  n is a positive integer input variable set to the number
#   of variables.
#
# x is an Abstractarray of length n. On input it must contain the
#   base point for the line search. On output it contains
#    x + stp * s.
#
# f is a variable. On input it must contain the value of f
#    at x. On output it contains the value of f at x + stp * s.
#
# g is an Abstractarray of length n. On input it must contain the
#    gradient of f at x. On output it contains the gradient
#    of f at x + stp * s.
#
# s is an input Abstractarray of length n which specifies the
#    search direction.
#
# stp is a nonnegative variable. On input stp contains an
#    initial estimate of a satisfactory step. On output
#    stp contains the final estimate.
#
#  f_tol and gtol are nonnegative input variables. Termination
#    occurs when the sufficient decrease condition and the
#    directional derivative condition are satisfied.
#
# x_tol is a nonnegative input variable. Termination occurs
#    when the relative width of the interval of uncertainty
#   is at most x_tol.
#
# alphamin and alphamax are nonnegative input variables which
#   specify lower and upper bounds for the step.
#
# maxfev is a positive integer input variable. Termination
#    occurs when the number of calls to df is at least
#    maxfev by the end of an iteration.
#
# info is an integer output variable set as follows:
#
#   info = 0  Improper input parameters.
#
#   info = 1  The sufficient decrease condition and the
#              directional derivative condition hold.
#
#   info = 2  Relative width of the interval of uncertainty
#            is at most x_tol.
#
#   info = 3  Number of calls to df has reached maxfev.
#
#   info = 4  The step is at the lower bound alphamin.
#
#   info = 5  The step is at the upper bound alphamax.
#
#   info = 6  Rounding errors prevent further progress.
#              There may not be a step which satisfies the
#              sufficient decrease and curvature conditions.
#              Tolerances may be too small.
#
#    nfev is an integer output variable set to the number of
#         calls to df.
#
#     Argonne National Laboratory. MINPACK Project. June 1983
#     Jorge J. More', David J. Thuente
#
#     **********

# Returns x, f, stp, info, nfev
# TODO: Decide whether to update x, f, g and info
#       or just return step and nfev and let existing code do its job

"""
The line search implementation from:
  Moré, Jorge J., and David J. Thuente
    Line search algorithms with guaranteed sufficient decrease.
    ACM Transactions on Mathematical Software (TOMS) 20.3 (1994): 286-307.
"""
@with_kw struct MoreThuente{T}
    f_tol::T = 1e-4 # c_1 Wolfe sufficient decrease condition
    gtol::T = 0.9   # c_2 Wolfe curvature condition (Recommend 0.1 for GradientDescent)
    x_tol::T = 1e-8
    alphamin::T = 1e-16
    alphamax::T = 65536.0
    maxfev::Int = 100
end

function (ls::MoreThuente)(df::AbstractObjective, x::AbstractArray{T},
                           s::AbstractArray{T}, alpha::Real, x_new::AbstractArray{T},
                           ϕ_0, dϕ_0) where T
    ϕdϕ = make_ϕdϕ(df, x_new, x, s)
    ls(ϕdϕ, alpha, ϕ_0, dϕ_0)
end

(ls::MoreThuente)(ϕ, dϕ, ϕdϕ, alpha, ϕ_0, dϕ_0) = ls(ϕdϕ, alpha, ϕ_0, dϕ_0)

# TODO: Should we deprecate the interface that only uses the ϕdϕ argument?
function (ls::MoreThuente)(ϕdϕ,
                           alpha::T,
                           ϕ_0,
                           dϕ_0) where T
    @unpack f_tol, gtol, x_tol, alphamin, alphamax, maxfev = ls

    iterfinitemax = -log2(eps(T))
    info = 0
    info_cstep = 1 # Info from step

    zeroT = convert(T, 0)

    #
    # Check the input parameters for errors.
    #

    if  alpha <= zeroT || f_tol < zeroT || gtol < zeroT ||
        x_tol < zeroT || alphamin < zeroT || alphamax < alphamin || maxfev <= zeroT
        throw(LineSearchException("Invalid parameters to MoreThuente.", 0))
    end

    if dϕ_0 >= zeroT
        throw(LineSearchException("Search direction is not a direction of descent.", 0))
    end

    #
    # Initialize local variables.
    #

    bracketed = false
    stage1 = true
    nfev = 0
    finit = ϕ_0
    dgtest = f_tol * dϕ_0
    width = alphamax - alphamin
    width1 = 2 * width

    # Keep this across calls

    #
    # The variables stx, fx, dgx contain the values of the step,
    # function, and directional derivative at the best step.
    # The variables sty, fy, dgy contain the value of the step,
    # function, and derivative at the other endpoint of
    # the interval of uncertainty.
    # The variables alpha, f, dg contain the values of the step,
    # function, and derivative at the current step.
    #

    stx = zeroT
    fx = finit
    dgx = dϕ_0
    sty = zeroT
    fy = finit
    dgy = dϕ_0

    # START: Ensure that the initial step provides finite function values
    # This is not part of the original FORTRAN code
    if !isfinite(alpha)
        alpha = one(T)
    end
    stmin = stx
    stmax = alpha + 4 * (alpha - stx) # Why 4?
    alpha = max(alpha, alphamin)
    alpha = min(alpha, alphamax)

    f, dg = ϕdϕ(alpha)
    nfev += 1 # This includes calls to f() and g!()
    iterfinite = 0
    while (!isfinite(f) || !isfinite(dg)) && iterfinite < iterfinitemax
        iterfinite += 1
        alpha = alpha/2

        f, dg = ϕdϕ(alpha)
        nfev += 1 # This includes calls to f() and g!()

        # Make stmax = (3/2)*alpha < 2alpha in the first iteration below
        stx = (convert(T, 7)/8)*alpha
    end
    # END: Ensure that the initial step provides finite function values

    while true
        #
        # Set the minimum and maximum steps to correspond
        # to the present interval of uncertainty.
        #

        if bracketed
            stmin = min(stx, sty)
            stmax = max(stx, sty)
        else
            stmin = stx
            stmax = alpha + 4 * (alpha - stx) # Why 4?
        end

        #
        # Ensure stmin and stmax (used in cstep) don't violate alphamin and alphamax
        # Not part of original FORTRAN translation
        #
        stmin = max(alphamin,stmin)
        stmax = min(alphamax,stmax)

        #
        # Force the step to be within the bounds alphamax and alphamin
        #

        alpha = max(alpha, alphamin)
        alpha = min(alpha, alphamax)

        #
        # If an unusual termination is to occur then let
        # alpha be the lowest point obtained so far.
        #

        if (bracketed && (alpha <= stmin || alpha >= stmax)) ||
            nfev >= maxfev-1 || info_cstep == 0 ||
            (bracketed && stmax - stmin <= x_tol * stmax)
            alpha = stx
        end

        #
        # Evaluate the function and gradient at alpha
        # and compute the directional derivative.
        #
        f, dg = ϕdϕ(alpha)
        nfev += 1 # This includes calls to f() and g!()

        if isapprox(dg, 0, atol=eps(T)) # Should add atol value to MoreThuente
            return alpha, f
        end

        ftest1 = finit + alpha * dgtest

        #
        # Test for convergence.
        #

        # What does info_cstep stand for?

        if (bracketed && (alpha <= stmin || alpha >= stmax)) || info_cstep == 0
            info = 6
        end
        if alpha == alphamax && f <= ftest1 && dg <= dgtest
            info = 5
        end
        if alpha == alphamin && (f > ftest1 || dg >= dgtest)
            info = 4
        end
        if nfev >= maxfev
            info = 3
        end
        if bracketed && stmax - stmin <= x_tol * stmax
            info = 2
        end
        if f <= ftest1 && abs(dg) <= -gtol * dϕ_0
            info = 1
        end

        #
        # Check for termination.
        #

        if info != 0
            return alpha, f
        end

        #
        # In the first stage we seek a step for which the modified
        # function has a nonpositive value and nonnegative derivative.
        #

        if stage1 && f <= ftest1 && dg >= min(f_tol, gtol) * dϕ_0
            stage1 = false
        end

        #
        # A modified function is used to predict the step only if
        # we have not obtained a step for which the modified
        # function has a nonpositive function value and nonnegative
        # derivative, and if a lower function value has been
        # obtained but the decrease is not sufficient.
        #

        if stage1 && f <= fx && f > ftest1
            #
            # Define the modified function and derivative values.
            #
            fm = f - alpha * dgtest
            fxm = fx - stx * dgtest
            fym = fy - sty * dgtest
            dgm = dg - dgtest
            dgxm = dgx - dgtest
            dgym = dgy - dgtest
            #
            # Call cstep to update the interval of uncertainty
            # and to compute the new step.
            #
            stx, fxm, dgxm,
            sty, fym, dgym,
            alpha, fm, dgm,
            bracketed, info_cstep =
                cstep(stx, fxm, dgxm, sty, fym, dgym,
                      alpha, fm, dgm, bracketed, stmin, stmax)
            #
            # Reset the function and gradient values for f.
            #
            fx = fxm + stx * dgtest
            fy = fym + sty * dgtest
            dgx = dgxm + dgtest
            dgy = dgym + dgtest
        else
            #
            # Call cstep to update the interval of uncertainty
            # and to compute the new step.
            #
            stx, fx, dgx,
            sty, fy, dgy,
            alpha, f, dg,
            bracketed, info_cstep =
                cstep(stx, fx, dgx, sty, fy, dgy,
                      alpha, f, dg, bracketed, stmin, stmax)
        end

        #
        # Force a sufficient decrease in the size of the
        # interval of uncertainty.
        #

        if bracketed
            if abs(sty - stx) >= (convert(T, 2)/3) * width1
                alpha = stx + (sty - stx)/2
            end
            width1 = width
            width = abs(sty - stx)
        end
    end # while
end # function


# Translation of minpack subroutine cstep
# Dianne O'Leary   July 1991
#
# Subroutine cstep
#
# The purpose of cstep is to compute a safeguarded step for
# a linesearch and to update an interval of uncertainty for
# a minimizer of the function.
#
# The parameter stx contains the step with the least function
# value. The parameter stp contains the current step. It is
# assumed that the derivative at stx is negative in the
# direction of the step. If bracketed is set true then a
# minimizer has been bracketed in an interval of uncertainty
# with endpoints stx and sty.
#
# The subroutine statement is
#
# subroutine cstep(stx, fx, dgx,
#                  sty, fy, dgy,
#                  stp, f, dg,
#                  bracketed, alphamin, alphamax, info)
#
# where
#
# stx, fx, and dgx are variables which specify the step,
#   the function, and the derivative at the best step obtained
#   so far. The derivative must be negative in the direction
#   of the step, that is, dgx and stp-stx must have opposite
#   signs. On output these parameters are updated appropriately
#
# sty, fy, and dgy are variables which specify the step,
#   the function, and the derivative at the other endpoint of
#   the interval of uncertainty. On output these parameters are
#   updated appropriately
#
# stp, f, and dg are variables which specify the step,
#   the function, and the derivative at the current step.
#   If bracketed is set true then on input stp must be
#   between stx and sty. On output stp is set to the new step
#
# bracketed is a logical variable which specifies if a minimizer
#   has been bracketed. If the minimizer has not been bracketed
#   then on input bracketed must be set false. If the minimizer
#   is bracketed then on output bracketed is set true
#
# alphamin and alphamax are input variables which specify lower
#   and upper bounds for the step
#
# info is an integer output variable set as follows:
#   If info = 1,2,3,4,5, then the step has been computed
#   according to one of the five cases below. Otherwise
#   info = 0, and this indicates improper input parameters
#
# Argonne National Laboratory. MINPACK Project. June 1983
# Jorge J. More', David J. Thuente

function cstep(stx::Real, fx::Real, dgx::Real,
               sty::Real, fy::Real, dgy::Real,
               alpha::Real, f::Real, dg::Real,
               bracketed::Bool, alphamin::Real, alphamax::Real)

   T = promote_type(typeof(stx), typeof(fx), typeof(dgx), typeof(sty), typeof(fy), typeof(dgy), typeof(alpha), typeof(f), typeof(dg), typeof(alphamin), typeof(alphamax))
   zeroT = convert(T, 0)
   info = 0

   #
   # Check the input parameters for error
   #

   if (bracketed && (alpha <= min(stx, sty) || alpha >= max(stx, sty))) ||
     dgx * (alpha - stx) >= zeroT || alphamax < alphamin
       throw(ArgumentError("Minimizer not bracketed"))
   end

   #
   # Determine if the derivatives have opposite sign
   #

   sgnd = dg * (dgx / abs(dgx))

   #
   # First case. A higher function value.
   # The minimum is bracketed. If the cubic step is closer
   # to stx than the quadratic step, the cubic step is taken,
   # else the average of the cubic and quadratic steps is taken
   #

   if f > fx
      info = 1
      bound = true
      theta = 3 * (fx - f) / (alpha - stx) + dgx + dg
      # Use s to prevent overflow/underflow of theta^2 and dgx * dg
      s = max(abs(theta), abs(dgx), abs(dg))
      gamma = s * sqrt((theta / s)^2 - (dgx / s) * (dg / s))
      if alpha < stx
          gamma = -gamma
      end
      p = gamma - dgx + theta
      q = gamma - dgx + gamma + dg
      r = p / q
      alphac = stx + r * (alpha - stx)
      alphaq = stx + (dgx / ((fx - f) / (alpha - stx) + dgx)) / 2 * (alpha - stx)
      if abs(alphac - stx) < abs(alphaq - stx)
         alphaf = alphac
      else
         alphaf = (alphac + alphaq) / 2
      end
      bracketed = true

   #
   # Second case. A lower function value and derivatives of
   # opposite sign. The minimum is bracketed. If the cubic
   # step is closer to stx than the quadratic (secant) step,
   # the cubic step is taken, else the quadratic step is taken
   #

elseif sgnd < zeroT
      info = 2
      bound = false
      theta = 3 * (fx - f) / (alpha - stx) + dgx + dg
      # Use s to prevent overflow/underflow of theta^2 and dgx * dg
      s = max(abs(theta), abs(dgx), abs(dg))
      gamma = s * sqrt((theta / s)^2 - (dgx / s) * (dg / s))

      if alpha > stx
         gamma = -gamma
      end
      p = gamma - dg + theta
      q = gamma - dg + gamma + dgx
      r = p / q
      alphac = alpha + r * (stx - alpha)
      alphaq = alpha + (dg / (dg - dgx)) * (stx - alpha)
      if abs(alphac - alpha) > abs(alphaq - alpha)
         alphaf = alphac
      else
         alphaf = alphaq
      end
      bracketed = true

   #
   # Third case. A lower function value, derivatives of the
   # same sign, and the magnitude of the derivative decreases.
   # The cubic step is only used if the cubic tends to infinity
   # in the direction of the step or if the minimum of the cubic
   # is beyond alpha. Otherwise the cubic step is defined to be
   # either alphamin or alphamax. The quadratic (secant) step is also
   # computed and if the minimum is bracketed then the the step
   # closest to stx is taken, else the step farthest away is taken
   #

   elseif abs(dg) < abs(dgx)
      info = 3
      bound = true
      theta = 3 * (fx - f) / (alpha - stx) + dgx + dg
      # Use s to prevent overflow/underflow of theta^2 and dgx * dg
      s = max(abs(theta), abs(dgx), abs(dg))
      #
      # The case gamma = 0 only arises if the cubic does not tend
      # to infinity in the direction of the step
      #
      # # Use NaNMath in case s == zero(s)
      gamma = s * sqrt(NaNMath.max(zero(s), (theta / s)^2 - (dgx / s) * (dg / s)))

      if alpha > stx
          gamma = -gamma
      end
      p = gamma - dg + theta
      q = gamma + dgx - dg + gamma
      r = p / q
      if r < zeroT && gamma != zeroT
         alphac = alpha + r * (stx - alpha)
     elseif alpha > stx
         alphac = alphamax
      else
         alphac = alphamin
      end
      alphaq = alpha + (dg / (dg - dgx)) * (stx - alpha)
      if bracketed
         if abs(alpha - alphac) < abs(alpha - alphaq)
            alphaf = alphac
         else
            alphaf = alphaq
         end
      else
         if abs(alpha - alphac) > abs(alpha - alphaq)
            alphaf = alphac
         else
            alphaf = alphaq
         end
      end

   #
   # Fourth case. A lower function value, derivatives of the
   # same sign, and the magnitude of the derivative does
   # not decrease. If the minimum is not bracketed, the step
   # is either alphamin or alphamax, else the cubic step is taken
   #

   else
      info = 4
      bound = false
      if bracketed
         theta = 3 * (f - fy) / (sty - alpha) + dgy + dg
         # Use s to prevent overflow/underflow of theta^2 and dgy * dg
         s = max(abs(theta), abs(dgy), abs(dg))
         gamma = s * sqrt((theta / s)^2 - (dgy / s) * (dg / s))

         if alpha > sty
             gamma = -gamma
         end
         p = gamma - dg + theta
         q = gamma - dg + gamma + dgy
         r = p / q
         alphac = alpha + r * (sty - alpha)
         alphaf = alphac
     elseif alpha > stx
         alphaf = alphamax
      else
         alphaf = alphamin
      end
   end

   #
   # Update the interval of uncertainty. This update does not
   # depend on the new step or the case analysis above
   #

   if f > fx
      sty = alpha
      fy = f
      dgy = dg
   else
      if sgnd < zeroT
         sty = stx
         fy = fx
         dgy = dgx
      end
      stx = alpha
      fx = f
      dgx = dg
   end

   #
   # Compute the new step and safeguard it
   #

   alphaf = min(alphamax, alphaf)
   alphaf = max(alphamin, alphaf)
   alpha = alphaf
   if bracketed && bound
      if sty > stx
         alpha = min(stx + (convert(T, 2)/3) * (sty - stx), alpha)
      else
         alpha = max(stx + (convert(T, 2)/3) * (sty - stx), alpha)
      end
   end

   return stx, fx, dgx, sty, fy, dgy, alpha, f, dg, bracketed, info
end
