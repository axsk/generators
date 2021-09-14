module EAMC

include("main.jl")
include("galerkin.jl")
include("perronfrobenius.jl")
include("commitor.jl")
include("processes/processes.jl")
#include("plots.jl")
include("utils.jl")
include("gillespie.jl")

export sqra, galerkin, commitor

end