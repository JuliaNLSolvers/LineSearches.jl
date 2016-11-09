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
#        subroutine cvsrch(df,n,x,f,g,s,stp,f_tol,gtol,x_tol,
#                          stpmin,stpmax,maxfev,info,nfev,wa)
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
# x is an array of length n. On input it must contain the
#   base point for the line search. On output it contains
#    x + stp * s.
#
# f is a variable. On input it must contain the value of f
#    at x. On output it contains the value of f at x + stp * s.
#
# g is an array of length n. On input it must contain the
#    gradient of f at x. On output it contains the gradient
#    of f at x + stp * s.
#
# s is an input array of length n which specifies the
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
# stpmin and stpmax are nonnegative input variables which
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
#   info = 4  The step is at the lower bound stpmin.
#
#   info = 5  The step is at the upper bound stpmax.
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

# Returns x, f, g, stp, info, nfev
# TODO: Decide whether to update x, f, g and info
#       or just return step and nfev and let existing code do its job

function morethuente!{T}(df,
                         x::Vector,
                         s::Vector,
                         new_x::Vector,
                         g::Vector,
                         lsr::LineSearchResults{T},
                         c::Real,
                         mayterminate::Bool;
                         n::Integer = length(x),
                         stp::Real = 1.0,
                         f_tol::Real = 1e-4,
                         gtol::Real = 0.9,
                         x_tol::Real = 1e-8,
                         stpmin::Real = 1e-16,
                         stpmax::Real = 65536.0,
                         maxfev::Integer = 100)

    info = 0
    info_cstep = 1 # Info from step

    #
    # Check the input parameters for errors.
    #

    if n <= 0 || stp <= 0.0 || f_tol < 0.0 || gtol < 0.0 ||
        x_tol < 0.0 || stpmin < 0.0 || stpmax < stpmin || maxfev <= 0
        throw(ArgumentError("Invalid parameters to morethuente"))
    end

    # Count function and gradient calls
    f_calls = 0
    g_calls = 0

    # read finit and slope from LineSearchResults
    f = lsr.value[end]
    dginit = lsr.slope[end] # dot(g,s)
    if dginit >= 0.0
        throw(ArgumentError("Search direction is not a direction of descent"))
    end

    #
    # Initialize local variables.
    #

    bracketed = false
    stage1 = true
    nfev = 0
    finit = f
    dgtest = f_tol * dginit
    width = stpmax - stpmin
    width1 = 2 * width
    copy!(new_x, x)
    # Keep this across calls
    # Replace with new_x?

    #
    # The variables stx, fx, dgx contain the values of the step,
    # function, and directional derivative at the best step.
    # The variables sty, fy, dgy contain the value of the step,
    # function, and derivative at the other endpoint of
    # the interval of uncertainty.
    # The variables stp, f, dg contain the values of the step,
    # function, and derivative at the current step.
    #

    stx = 0.0
    fx = finit
    dgx = dginit
    sty = 0.0
    fy = finit
    dgy = dginit

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
            stmax = stp + 4 * (stp - stx) # Why 4?
        end

        #
        # Force the step to be within the bounds stpmax and stpmin
        #

        stp = max(stp, stpmin)
        stp = min(stp, stpmax)

        #
        # If an unusual termination is to occur then let
        # stp be the lowest point obtained so far.
        #

        if (bracketed && (stp <= stmin || stp >= stmax)) ||
            nfev >= maxfev-1 || info_cstep == 0 ||
            (bracketed && stmax - stmin <= x_tol * stmax)
            stp = stx
        end

        #
        # Evaluate the function and gradient at stp
        # and compute the directional derivative.
        #

        for i in 1:n
            new_x[i] = x[i] + stp * s[i] # TODO: Use x_new here
        end
        f = df.fg!(new_x, g)
        f_calls += 1
        g_calls += 1
        nfev += 1 # This includes calls to f() and g!()
        dg = vecdot(g, s)
        ftest1 = finit + stp * dgtest

        #
        # Test for convergence.
        #

        # What does info_cstep stand for?

        if (bracketed && (stp <= stmin || stp >= stmax)) || info_cstep == 0
            info = 6
        end
        if stp == stpmax && f <= ftest1 && dg <= dgtest
            info = 5
        end
        if stp == stpmin && (f > ftest1 || dg >= dgtest)
            info = 4
        end
        if nfev >= maxfev
            info = 3
        end
        if bracketed && stmax - stmin <= x_tol * stmax
            info = 2
        end
        if f <= ftest1 && abs(dg) <= -gtol * dginit
            info = 1
        end

        #
        # Check for termination.
        #

        if info != 0
            return stp, f_calls, g_calls
        end

        #
        # In the first stage we seek a step for which the modified
        # function has a nonpositive value and nonnegative derivative.
        #

        if stage1 && f <= ftest1 && dg >= min(f_tol, gtol) * dginit
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
            fm = f - stp * dgtest
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
            stp, fm, dgm,
            bracketed, info_cstep =
                cstep(stx, fxm, dgxm, sty, fym, dgym,
                      stp, fm, dgm, bracketed, stmin, stmax)
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
            stp, f, dg,
            bracketed, info_cstep =
                cstep(stx, fx, dgx, sty, fy, dgy,
                      stp, f, dg, bracketed, stmin, stmax)
        end

        #
        # Force a sufficient decrease in the size of the
        # interval of uncertainty.
        #

        if bracketed
            if abs(sty - stx) >= 0.66 * width1
                stp = stx + 0.5 * (sty - stx)
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
#                  stp, f, dp,
#                  bracketed, stpmin, stpmax, info)
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
# stp, f, and dp are variables which specify the step,
#   the function, and the derivative at the current step.
#   If bracketed is set true then on input stp must be
#   between stx and sty. On output stp is set to the new step
#
# bracketed is a logical variable which specifies if a minimizer
#   has been bracketed. If the minimizer has not been bracketed
#   then on input bracketed must be set false. If the minimizer
#   is bracketed then on output bracketed is set true
#
# stpmin and stpmax are input variables which specify lower
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
               stp::Real, f::Real, dg::Real,
               bracketed::Bool, stpmin::Real, stpmax::Real)

   info = 0

   #
   # Check the input parameters for error
   #

   if (bracketed && (stp <= min(stx, sty) || stp >= max(stx, sty))) ||
        dgx * (stp - stx) >= 0.0 || stpmax < stpmin
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
      theta = 3 * (fx - f) / (stp - stx) + dgx + dg
      s = max(theta, dgx, dg)
      gamma = s * sqrt((theta / s)^2 - (dgx / s) * (dg / s))
      if stp < stx
          gamma = -gamma
      end
      p = gamma - dgx + theta
      q = gamma - dgx + gamma + dg
      r = p / q
      stpc = stx + r * (stp - stx)
      stpq = stx + ((dgx / ((fx - f) / (stp - stx) + dgx)) / 2) * (stp - stx)
      if abs(stpc - stx) < abs(stpq - stx)
         stpf = stpc
      else
         stpf = stpc + (stpq - stpc) / 2
      end
      bracketed = true

   #
   # Second case. A lower function value and derivatives of
   # opposite sign. The minimum is bracketed. If the cubic
   # step is closer to stx than the quadratic (secant) step,
   # the cubic step is taken, else the quadratic step is taken
   #

   elseif sgnd < 0.0
      info = 2
      bound = false
      theta = 3 * (fx - f) / (stp - stx) + dgx + dg
      s = max(theta, dgx, dg)
      gamma = s * sqrt((theta / s)^2 - (dgx / s) * (dg / s))
      if stp > stx
         gamma = -gamma
      end
      p = gamma - dg + theta
      q = gamma - dg + gamma + dgx
      r = p / q
      stpc = stp + r * (stx - stp)
      stpq = stp + (dg / (dg - dgx)) * (stx - stp)
      if abs(stpc - stp) > abs(stpq - stp)
         stpf = stpc
      else
         stpf = stpq
      end
      bracketed = true

   #
   # Third case. A lower function value, derivatives of the
   # same sign, and the magnitude of the derivative decreases.
   # The cubic step is only used if the cubic tends to infinity
   # in the direction of the step or if the minimum of the cubic
   # is beyond stp. Otherwise the cubic step is defined to be
   # either stpmin or stpmax. The quadratic (secant) step is also
   # computed and if the minimum is bracketed then the the step
   # closest to stx is taken, else the step farthest away is taken
   #

   elseif abs(dg) < abs(dgx)
      info = 3
      bound = true
      theta = 3 * (fx - f) / (stp - stx) + dgx + dg
      s = max(theta, dgx, dg)
      #
      # The case gamma = 0 only arises if the cubic does not tend
      # to infinity in the direction of the step
      #
      gamma = s * sqrt(max(0.0, (theta / s)^2 - (dgx / s) * (dg / s)))
      if stp > stx
          gamma = -gamma
      end
      p = gamma - dg + theta
      q = gamma + dgx - dg + gamma
      r = p / q
      if r < 0.0 && gamma != 0.0
         stpc = stp + r * (stx - stp)
      elseif stp > stx
         stpc = stpmax
      else
         stpc = stpmin
      end
      stpq = stp + (dg / (dg - dgx)) * (stx - stp)
      if bracketed
         if abs(stp - stpc) < abs(stp - stpq)
            stpf = stpc
         else
            stpf = stpq
         end
      else
         if abs(stp - stpc) > abs(stp - stpq)
            stpf = stpc
         else
            stpf = stpq
         end
      end

   #
   # Fourth case. A lower function value, derivatives of the
   # same sign, and the magnitude of the derivative does
   # not decrease. If the minimum is not bracketed, the step
   # is either stpmin or stpmax, else the cubic step is taken
   #

   else
      info = 4
      bound = false
      if bracketed
         theta = 3 * (f - fy) / (sty - stp) + dgy + dg
         s = max(theta, dgy, dg)
         gamma = s * sqrt((theta / s)^2 - (dgy / s) * (dg / s))
         if stp > sty
             gamma = -gamma
         end
         p = gamma - dg + theta
         q = gamma - dg + gamma + dgy
         r = p / q
         stpc = stp + r * (sty - stp)
         stpf = stpc
      elseif stp > stx
         stpf = stpmax
      else
         stpf = stpmin
      end
   end

   #
   # Update the interval of uncertainty. This update does not
   # depend on the new step or the case analysis above
   #

   if f > fx
      sty = stp
      fy = f
      dgy = dg
   else
      if sgnd < 0.0
         sty = stx
         fy = fx
         dgy = dgx
      end
      stx = stp
      fx = f
      dgx = dg
   end

   #
   # Compute the new step and safeguard it
   #

   stpf = min(stpmax, stpf)
   stpf = max(stpmin, stpf)
   stp = stpf
   if bracketed && bound
      if sty > stx
         stp = min(stx + 0.66 * (sty - stx), stp)
      else
         stp = max(stx + 0.66 * (sty - stx), stp)
      end
   end

   return stx, fx, dgx, sty, fy, dgy, stp, f, dg, bracketed, info
end
