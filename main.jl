include("processes/processes.jl")
include("embeddedspacetime.jl")
include("pccap.jl")

function cluster(process, timegrid; adj=true)
    G = augmentedembeddedmatrix(process, timegrid)
    heatmap(G) |> display
    G = (adj ? G' : G) |> collect
    e,v = sortschur(schur(G))
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

begin
    G,e,v = cluster(BarrierSwitch(), [0,1], adj=true)
    plot(heatmap(v), heatmap(e), layout=grid(2,1))
end

function massageev(v,e)

end
