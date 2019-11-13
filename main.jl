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

if false 
    G,e,v = cluster(BarrierSwitch(), [0,1], adj=true)
    plot(heatmap(v), heatmap(e), layout=grid(2,1))
end

gauss(x, mu, sigma) = 1 / (sigma * sqrt(2*pi)) * exp( - 1/2 * ((x-mu)/sigma)^2)

PotentialShiftingBarrier() = (x,t) -> gauss(x, 0.2 + t/2, 0.5)

function experiment(V, xs, ts, comittor, flux=1)
    p = OverdampedLangevin(V, flux, xs, false)
    g = galerkin(p, ts)

    term = zeros(length(xs))
    term[ceil(length(xs)/2):end] .= 1
    
    if commitor != nothing
        c = commitor(g, term)
    end
    
    Base.@locals
end

