include("processes/processes.jl")

include("embeddedspacetime.jl")
function test(; beta=1, flux=1, tmax=1)

    process = DoubleWell(beta, flux, 20)
    timegrid = range(0, stop=tmax, length=3)
    G = augmentedembeddedmatrix(process, timegrid)
    heatmap(G) |> display
    e,v = sortschur(schur(G' |> collect))
    @show e
    plot(v[:,1:6])

end

test(flux=3)
