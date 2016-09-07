module LineSearches

export AbstractDifferentiableFunction, LineSearchResults

export clear!, alphatry, alphainit

export hz_linesearch!, backtracking_linesearch!, interpolating_linesearch!,
    mt_linesearch!

include("types.jl")
include("api.jl")

# Line Search Methods
include("backtracking_linesearch.jl")
include("interpolating_linesearch.jl")
include("mt_cstep.jl")
include("mt_linesearch.jl")
include("hz_linesearch.jl")


end # module
