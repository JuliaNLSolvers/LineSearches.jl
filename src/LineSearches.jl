isdefined(Base, :__precompile__) && __precompile__()

module LineSearches

export LineSearchResults

export clear!, alphatry, alphainit

export hagerzhang!, backtracking!, strongwolfe!,
    morethuente!, interpbacktrack!

include("types.jl")
include("api.jl")

# Line Search Methods
include("backtracking_linesearch.jl")
include("interpolating_linesearch.jl")
include("mt_cstep.jl")
include("mt_linesearch.jl")
include("hz_linesearch.jl")


end # module
