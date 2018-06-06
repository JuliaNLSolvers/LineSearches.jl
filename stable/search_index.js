var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "DocTestSetup = :(using LineSearches)"
},

{
    "location": "index.html#LineSearches.jl-1",
    "page": "Home",
    "title": "LineSearches.jl",
    "category": "section",
    "text": "A line search toolbox written in Julia."
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "LineSearches provides a collection of line search routines for optimization and nonlinear solvers.  The package can be used on its own, but it also provides extra supporting functionality for Optim.jl and NLsolve.jl."
},

{
    "location": "index.html#Available-line-search-algorithms-1",
    "page": "Home",
    "title": "Available line search algorithms",
    "category": "section",
    "text": "HagerZhang (Taken from the Conjugate Gradient implementation by Hager and Zhang, 2006)\nMoreThuente (From the algorithm in More and Thuente, 1994)\nBackTracking (Described in Nocedal and Wright, 2006)\nStrongWolfe (Nocedal and Wright)\nStatic (Takes the proposed initial step length.)"
},

{
    "location": "index.html#Available-initial-step-length-procedures-1",
    "page": "Home",
    "title": "Available initial step length procedures",
    "category": "section",
    "text": "The package provides some procedures to calculate the initial step length that is passed to the line search algorithm, currently specialized to be used with Optim and NLsolve.InitialPrevious (Use the step length from the previous optimization iteration)\nInitialStatic (Use the same initial step length each time)\nInitialHagerZhang (Taken from Hager and Zhang, 2006)\nInitialQuadratic (Propose initial step length based on a quadratic interpolation)\nInitialConstantChange (Propose initial step length assuming constant change in step length)"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "To install, simply run the following in the Julia REPL:Pkg.add(\"LineSearches\")and then runusing LineSearchesto load the package."
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "W. W. Hager and H. Zhang (2006) \"Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed descent.\" ACM Transactions on Mathematical Software 32: 113-137.\nMoré, Jorge J., and David J. Thuente. \"Line search algorithms with guaranteed sufficient decrease.\" ACM Transactions on Mathematical Software (TOMS) 20.3 (1994): 286-307.\nNocedal, Jorge, and Stephen Wright. \"Numerical optimization.\" Springer Science & Business Media, 2006."
},

{
    "location": "examples/generated/customoptimizer.html#",
    "page": "Using LineSearches without Optim/NLsolve",
    "title": "Using LineSearches without Optim/NLsolve",
    "category": "page",
    "text": "EditURL = \"https://github.com/JuliaNLSolvers/LineSearches.jl/blob/master/docs/src/examples/customoptimizer.jl\""
},

{
    "location": "examples/generated/customoptimizer.html#Using-LineSearches-without-Optim/NLsolve-1",
    "page": "Using LineSearches without Optim/NLsolve",
    "title": "Using LineSearches without Optim/NLsolve",
    "category": "section",
    "text": "tip: Tip\nThis example is also available as a Jupyter notebook: customoptimizer.ipynbThis tutorial shows you how to use the line search algorithms in LineSearches for your own optimization algorithm that is not part of Optim or NLsolve.Say we have written a gradient descent optimization algorithm but would like to experiment with different line search algorithms. The algorithm is implemented as follows.function gdoptimize(f, g!, fg!, x0::AbstractArray{T}, linesearch,\n                    maxiter::Int = 10000,\n                    g_rtol::T = sqrt(eps(T)), g_atol::T = eps(T)) where T <: Number\n    x = copy(x0)\n    gvec = similar(x)\n    g!(gvec, x)\n    fx = f(x)\n\n    gnorm = norm(gvec)\n    gtol = max(g_rtol*gnorm, g_atol)\n\n    # Univariate line search functions\n    ϕ(α) = f(x .+ α.*s)\n    function dϕ(α)\n        g!(gvec, x .+ α.*s)\n        return vecdot(gvec, s)\n    end\n    function ϕdϕ(α)\n        phi = fg!(gvec, x .+ α.*s)\n        dphi = vecdot(gvec, s)\n        return (phi, dphi)\n    end\n\n    s = similar(gvec) # Step direction\n\n    iter = 0\n    while iter < maxiter && gnorm > gtol\n        iter += 1\n        s .= -gvec\n\n        dϕ_0 = dot(s, gvec)\n        α, fx = linesearch(ϕ, dϕ, ϕdϕ, 1.0, fx, dϕ_0)\n\n        @. x = x + α*s\n        g!(gvec, x)\n        gnorm = norm(gvec)\n    end\n\n    return (fx, x, iter)\nendNote that there are many optimization and line search algorithms that allow the user to evaluate both the objective and the gradient at the same time, for computational efficiency reasons. We have included this functionality in the algorithm as the input function fg!, and even if the Gradient Descent algorithm does not use it explicitly, many of the LineSearches algorithms do.The Gradient Descent gdoptimize method selects a descent direction and calls the line search algorithm  linesearch which returns the step length α and the objective value fx = f(x + α*s).The functions ϕ and dϕ represent a univariate objective and its derivative, which is used by the line search algorithms. To utilize the fg! function call in the optimizer, some of the line searches require a function ϕdϕ which returns the univariate objective and the derivative at the same time."
},

{
    "location": "examples/generated/customoptimizer.html#Optimizing-Rosenbrock-1",
    "page": "Using LineSearches without Optim/NLsolve",
    "title": "Optimizing Rosenbrock",
    "category": "section",
    "text": "Here is an example to show how we can combine gdoptimize and LineSearches to minimize the Rosenbrock function, which is defined byf(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2\n\nfunction g!(gvec, x)\n    gvec[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]\n    gvec[2] = 200.0 * (x[2] - x[1]^2)\n    gvec\nend\n\nfunction fg!(gvec, x)\n    g!(gvec, x)\n    f(x)\nendWe can now use gdoptimize with BackTracking to optimize the Rosenbrock function from a given initial condition x0.x0 = [-1., 1.0]\n\nusing LineSearches\nls = BackTracking(order=3)\nfx_bt3, x_bt3, iter_bt3 = gdoptimize(f, g!, fg!, x0, ls)Interestingly, the StrongWolfe line search converges in one iteration, whilst all the other algorithms take thousands of iterations. This is just luck due to the particular choice of initial conditionls = StrongWolfe()\nfx_sw, x_sw, iter_sw = gdoptimize(f, g!, fg!, x0, ls)"
},

{
    "location": "examples/generated/customoptimizer.html#customoptimizer-plain-program-1",
    "page": "Using LineSearches without Optim/NLsolve",
    "title": "Plain Program",
    "category": "section",
    "text": "Below follows a version of the program without any comments. The file is also available here: customoptimizer.jlfunction gdoptimize(f, g!, fg!, x0::AbstractArray{T}, linesearch,\n                    maxiter::Int = 10000,\n                    g_rtol::T = sqrt(eps(T)), g_atol::T = eps(T)) where T <: Number\n    x = copy(x0)\n    gvec = similar(x)\n    g!(gvec, x)\n    fx = f(x)\n\n    gnorm = norm(gvec)\n    gtol = max(g_rtol*gnorm, g_atol)\n\n    # Univariate line search functions\n    ϕ(α) = f(x .+ α.*s)\n    function dϕ(α)\n        g!(gvec, x .+ α.*s)\n        return vecdot(gvec, s)\n    end\n    function ϕdϕ(α)\n        phi = fg!(gvec, x .+ α.*s)\n        dphi = vecdot(gvec, s)\n        return (phi, dphi)\n    end\n\n    s = similar(gvec) # Step direction\n\n    iter = 0\n    while iter < maxiter && gnorm > gtol\n        iter += 1\n        s .= -gvec\n\n        dϕ_0 = dot(s, gvec)\n        α, fx = linesearch(ϕ, dϕ, ϕdϕ, 1.0, fx, dϕ_0)\n\n        @. x = x + α*s\n        g!(gvec, x)\n        gnorm = norm(gvec)\n    end\n\n    return (fx, x, iter)\nend\n\nf(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2\n\nfunction g!(gvec, x)\n    gvec[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]\n    gvec[2] = 200.0 * (x[2] - x[1]^2)\n    gvec\nend\n\nfunction fg!(gvec, x)\n    g!(gvec, x)\n    f(x)\nend\n\nx0 = [-1., 1.0]\n\nusing LineSearches\nls = BackTracking(order=3)\nfx_bt3, x_bt3, iter_bt3 = gdoptimize(f, g!, fg!, x0, ls)\n\nls = StrongWolfe()\nfx_sw, x_sw, iter_sw = gdoptimize(f, g!, fg!, x0, ls)\n\n# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jlThis page was generated using Literate.jl."
},

{
    "location": "examples/generated/optim_linesearch.html#",
    "page": "Optim line search",
    "title": "Optim line search",
    "category": "page",
    "text": ""
},

{
    "location": "examples/generated/optim_linesearch.html#Optim-line-search-1",
    "page": "Optim line search",
    "title": "Optim line search",
    "category": "section",
    "text": "tip: Tip\nThis example is also available as a Jupyter notebook: optim_linesearch.ipynbThis example shows how to use LineSearches with Optim.  We solve the Rosenbrock problem with two different line search algorithms.First, run Newton with the default line search algorithm:using Optim, LineSearches\nimport OptimTestProblems.MultivariateProblems\nUP = MultivariateProblems.UnconstrainedProblems\nprob = UP.examples[\"Rosenbrock\"]\n\nalgo_hz = Newton(linesearch = HagerZhang())\nres_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)Now we can try Newton with the cubic backtracking line search, which reduced the number of objective and gradient calls.algo_bt3 = Newton(linesearch = BackTracking(order=3))\nres_bt3 = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_bt3)This page was generated using Literate.jl."
},

{
    "location": "examples/generated/optim_initialstep.html#",
    "page": "Optim initial step length guess",
    "title": "Optim initial step length guess",
    "category": "page",
    "text": ""
},

{
    "location": "examples/generated/optim_initialstep.html#Optim-initial-step-length-guess-1",
    "page": "Optim initial step length guess",
    "title": "Optim initial step length guess",
    "category": "section",
    "text": "tip: Tip\nThis example is also available as a Jupyter notebook: optim_initialstep.ipynbThis example shows how to use the initial step length procedures with Optim.  We solve the Rosenbrock problem with two different procedures.First, run Newton with the (default) initial guess and line search procedures.using Optim, LineSearches\nimport OptimTestProblems.MultivariateProblems\nUP = MultivariateProblems.UnconstrainedProblems\nprob = UP.examples[\"Rosenbrock\"]\n\nalgo_st = Newton(alphaguess = InitialStatic(), linesearch = HagerZhang())\nres_st = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_st)We can now try with the initial step length guess from Hager and Zhang.algo_hz = Newton(alphaguess = InitialHagerZhang(α0=1.0), linesearch = HagerZhang())\nres_hz = Optim.optimize(prob.f, prob.g!, prob.h!, prob.initial_x, method=algo_hz)From the result we see that this has reduced the number of function and gradient calls, but increased the number of iterations.This page was generated using Literate.jl."
},

{
    "location": "reference/linesearch.html#",
    "page": "Line search routines",
    "title": "Line search routines",
    "category": "page",
    "text": "DocTestSetup = :(using LineSearches)"
},

{
    "location": "reference/linesearch.html#LineSearches.BackTracking",
    "page": "Line search routines",
    "title": "LineSearches.BackTracking",
    "category": "type",
    "text": "BackTracking specifies a backtracking line-search that uses a quadratic or cubic interpolant to determine the reduction in step-size. E.g., if f(α) > f(0) + c₁ α f\'(0), then the quadratic interpolant of f(0), f\'(0), f(α) has a minimiser α\' in the open interval (0, α). More strongly, there exists a factor ρ = ρ(c₁) such that α\' ≦ ρ α.\n\nThis is a modification of the algorithm described in Nocedal Wright (2nd ed), Sec. 3.5.\n\n\n\n"
},

{
    "location": "reference/linesearch.html#LineSearches.HagerZhang",
    "page": "Line search routines",
    "title": "LineSearches.HagerZhang",
    "category": "type",
    "text": "Conjugate gradient line search implementation from:   W. W. Hager and H. Zhang (2006) Algorithm 851: CG_DESCENT, a     conjugate gradient method with guaranteed descent. ACM     Transactions on Mathematical Software 32: 113–137.\n\n\n\n"
},

{
    "location": "reference/linesearch.html#LineSearches.MoreThuente",
    "page": "Line search routines",
    "title": "LineSearches.MoreThuente",
    "category": "type",
    "text": "The line search implementation from:   Moré, Jorge J., and David J. Thuente     Line search algorithms with guaranteed sufficient decrease.     ACM Transactions on Mathematical Software (TOMS) 20.3 (1994): 286-307.\n\n\n\n"
},

{
    "location": "reference/linesearch.html#LineSearches.Static",
    "page": "Line search routines",
    "title": "LineSearches.Static",
    "category": "function",
    "text": "Static: defines a static linesearch which returns the initial step length.\n\nStatic is intended for methods with well-scaled updates; i.e. Newton, on well-behaved problems.\n\n\n\n"
},

{
    "location": "reference/linesearch.html#LineSearches.StrongWolfe",
    "page": "Line search routines",
    "title": "LineSearches.StrongWolfe",
    "category": "type",
    "text": "StrongWolfe: This linesearch algorithm guarantees that the step length satisfies the (strong) Wolfe conditions. See Nocedal and Wright - Algorithms 3.5 and 3.6\n\nThis algorithm is mostly of theoretical interest, users should most likely use MoreThuente, HagerZhang or BackTracking.\n\nParameters:  (and defaults)\n\nc_1 = 1e-4: Armijo condition\nc_2 = 0.9 : second (strong) Wolfe condition\nρ = 2.0 : bracket growth\n\n\n\n"
},

{
    "location": "reference/linesearch.html#Line-search-routines-1",
    "page": "Line search routines",
    "title": "Line search routines",
    "category": "section",
    "text": "BackTracking\nHagerZhang\nMoreThuente\nStatic\nStrongWolfe"
},

{
    "location": "reference/initialstep.html#",
    "page": "Initial step length guess",
    "title": "Initial step length guess",
    "category": "page",
    "text": "DocTestSetup = :(using LineSearches)"
},

{
    "location": "reference/initialstep.html#LineSearches.InitialStatic",
    "page": "Initial step length guess",
    "title": "LineSearches.InitialStatic",
    "category": "type",
    "text": "Provide static initial step length.\n\nKeyword alpha corresponds to static step length, default is 1.0. If keyword scaled = true, then the initial step length is scaled with the l_2 norm of the step direction.\n\n\n\n"
},

{
    "location": "reference/initialstep.html#LineSearches.InitialPrevious",
    "page": "Initial step length guess",
    "title": "LineSearches.InitialPrevious",
    "category": "type",
    "text": "Use previous step length as initial guess, within the bounds [alphamin, alphamax]\n\nIf state.alpha is NaN, then return fallback value is.alpha\n\n\n\n"
},

{
    "location": "reference/initialstep.html#LineSearches.InitialQuadratic",
    "page": "Initial step length guess",
    "title": "LineSearches.InitialQuadratic",
    "category": "type",
    "text": "Quadratic interpolation for initial step length guess.\n\nThis is meant for methods that do not produce well-scaled search directions, such as Gradient Descent and (variations of) Conjugate Gradient methods. See the discussion around Nocedal and Wright, 2nd ed, (3.60).\n\nThis procedure have several arguments, with the following defaults.\n\nα0       = 1.0.         The initial step size at the first iteration.\nαmin     = 1e-12.       The minimum initial step size. (Default arbitrary).\nαmax     = 1.0.         The maximum initial step size.\nρ        = 0.25.        Maximum decrease from previous iteration, αinit ≥ α_{k-1}. (Default arbitrary).\nsnap2one = (0.75, Inf). Set all values within this (closed) interval to 1.0. (Default arbitrary).\n\nIf αmax ≠ 1.0, then you should consider to ensure that snap2one[2] < αmax.\n\n\n\n"
},

{
    "location": "reference/initialstep.html#LineSearches.InitialConstantChange",
    "page": "Initial step length guess",
    "title": "LineSearches.InitialConstantChange",
    "category": "type",
    "text": "Constant first-order change approximation to determine initial step length.\n\n** This requires that the optimization algorithm stores dphi0 from the previous iteration ** (dphi0_previous = real(vecdot(∇f_{k-1}, s_{k-1})), where s is the step direction.\n\nThis is meant for methods that do not produce well-scaled search directions, such as Gradient Descent and (variations of) Conjugate Gradient methods. See the discussion in Nocedal and Wright, 2nd ed, p. 59 on \"Initial Step Length\"\n\nThis procedure have several arguments, with the following defaults.\n\nα0       = 1.0.         The initial step size at the first iteration.\nαmin     = 1e-12.       The minimum initial step size. (Default arbitrary).\nαmax     = 1.0.         The maximum initial step size.\nρ        = 0.25.        Maximum decrease from previous iteration, αinit ≥ α_{k-1}. (Default arbitrary).\nsnap2one = (0.75, Inf). Set all values within this (closed) interval to 1.0. (Default arbitrary).\n\nIf αmax ≠ 1.0, then you should consider to ensure that snap2one[2] < αmax.\n\n\n\n"
},

{
    "location": "reference/initialstep.html#LineSearches.InitialHagerZhang",
    "page": "Initial step length guess",
    "title": "LineSearches.InitialHagerZhang",
    "category": "type",
    "text": "Initial step size algorithm from   W. W. Hager and H. Zhang (2006) Algorithm 851: CG_DESCENT, a     conjugate gradient method with guaranteed descent. ACM     Transactions on Mathematical Software 32: 113–137.\n\nIf α0 is NaN, then procedure I0 is called at the first iteration, otherwise, we select according to procedure I1-2, with starting value α0.\n\n\n\n"
},

{
    "location": "reference/initialstep.html#Initial-step-length-guess-1",
    "page": "Initial step length guess",
    "title": "Initial step length guess",
    "category": "section",
    "text": "Some of these routines are tightly integrated with Optim.InitialStatic\nInitialPrevious\nInitialQuadratic\nInitialConstantChange\nInitialHagerZhang"
},

]}
