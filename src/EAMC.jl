module EAMC

include("main.jl")
include("galerkin.jl")
include("processes/processes.jl")
include("plots.jl")
include("utils.jl")

export sqra, galerkin, commitor

end