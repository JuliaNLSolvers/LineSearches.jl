isdefined(Base, :__precompile__) && __precompile__()

module LineSearches

export LineSearchResults

export clear!, alphatry, alphainit

export hagerzhang!, backtracking!, strongwolfe!,
    morethuente!, interpbacktrack!

include("types.jl")
include("api.jl")

# Line Search Methods
include("backtracking.jl")
include("strongwolfe.jl")
include("mt_cstep.jl")
include("morethuente.jl")
include("hagerzhang.jl")


end # module
