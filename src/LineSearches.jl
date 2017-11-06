isdefined(Base, :__precompile__) && __precompile__()

module LineSearches

using Parameters, NaNMath

import NLSolversBase
import Base.clear!

export LineSearchResults, LineSearchException

export clear!, alphatry, alphainit

export BackTracking, HagerZhang, Static, MoreThuente, StrongWolfe

include("types.jl")
include("api.jl")

# Line Search Methods
include("backtracking.jl")
include("strongwolfe.jl")
include("morethuente.jl")
include("hagerzhang.jl")
include("static.jl")

include("deprecate.jl")

end # module
