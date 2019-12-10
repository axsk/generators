using Parameters


#  1d process with 4 cells divided by barriers (i.e. no exchange), where the central barrier is turned on at time `tswitch`
#  [ --- | --- (|) --- | ---]
@with_kw struct BarrierSwitch <: Process
    tswitch = [1] # switchting points, process starts closed. e.g. [0] for constant open
    cellresolution = 2
    rate = 1
end

length(p::BarrierSwitch) = 4 * p.cellresolution

# todo: add temperature/noise and flux-constant!
function generatormatrix(p::BarrierSwitch, t)
    @unpack_BarrierSwitch p

    c = cellresolution
    n = 4 * cellresolution
    G = zeros(n,n)

    # square root approximation

    for cell = 0:3
        for i = cell*cellresolution+1 : (cell+1)*cellresolution - 1
            G[i, i+1] = 1
            G[i+1, i] = 1
        end
    end

    if findfirst(x->x>=t, [tswitch...,Inf]) % 2 == 0
        i = 2 * cellresolution
        G[i, i+1] = 1
        G[i+1, i] = 1
    end

    for i = 1:n
        G[i,i] = -sum(G[i,:])
    end
    G * rate
end
