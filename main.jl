include("processes/processes.jl")
include("embeddedspacetime.jl")

function cluster(process, timegrid)
    G = augmentedembeddedmatrix(process, timegrid)
    heatmap(G) |> display
    e,v = sortschur(schur(G' |> collect))
    plot(v[:,1:6]) |> display
    G, e, v
end


function test1()
    process = DoubleWell(beta, flux, 20)
    timegrid = range(0, stop=tmax, length=3)
    cluster(process, timegrid)
end

function test2()
    process = BarrierSwitch()
    timegrid = [0,1]
    cluster(process, timegrid)
end
