isdefined(Base, :__precompile__) && __precompile__()

module LineSearches

using Parameters

import NLSolversBase
import Base.clear!

export LineSearchResults, LineSearchException

export clear!, alphatry, alphainit

export hagerzhang!, backtracking!, strongwolfe!,
    morethuente!, bt2!, bt3!, basic!,
    BackTracking, HagerZhang, Static, MoreThuente, StrongWolfe

include("types.jl")
include("api.jl")

# Line Search Methods
include("backtracking.jl")
include("strongwolfe.jl")
include("morethuente.jl")
include("hagerzhang.jl")
include("basic.jl")

include("deprecate.jl")

end # module
