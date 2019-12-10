abstract type Process end

import Base.length

include("langevin.jl")
include("barrierswitch.jl")
include("continuity.jl")
include("bickleyjet.jl")