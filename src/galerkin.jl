using Arpack, LinearAlgebra
using SparseArrays

export galerkin

#galerkin(Qs::Vector{Array}, dt) = galerkin(map(sparse, Qs), dt)
function galerkin(Qs::Vector{<:SparseMatrixCSC}, dt)
    n = size(Qs[1], 1)
    m = length(Qs)
    @assert length(dt) == m

    I = Int[]
    J = Int[]
    V = Float64[]

    qout = [collect(-diag(Qs[t])) for t in 1:m]
    s = [exp.(-dt[t]*qout[t]) for t in 1:m] # note that we still need the possible 0 rates here (c.f. below). in case of dt == Inf this can lead to 0*Inf = NaN
    for i=1:m
        replace!(qout[i], 0=>1.) # in case of 0 rates, replace with 1 to avoid division by 0. should not change any result since we multiply corresponding rows with zero rates later.
    end
    qt = [dropzeros((Qs[t] - Diagonal(Qs[t])) ./ qout[t]) for t in 1:m]

    for ti in 1:m
        for tj in ti:m
            fact = 1 ./ qout[ti] ./ dt[ti]
            if ti==tj
                fact2 = (s[ti] + dt[ti] * qout[ti] .- 1)
            elseif ti<tj
                fact2 = (1 .-s[ti]) .* (1 .-s[tj]) .* mycumprod(s, ti, tj)
            end
            res  = qt[tj] .* (fact .* fact2)
            II,JJ,VV = findnz(res)
            append!(I, II.+(ti-1)*n)
            append!(J, JJ.+(tj-1)*n)
            append!(V, VV)
        end
    end

    G = sparse(I,J,V, n*m, n*m)

    if dt[end] == Inf  # account for absorbing boundary
        G[(m-1)*n+1:end, (m-1)*n+1:end] = qt[m]
    end
    G
end

# cumulative product used in the holding probability between different time boxes
function mycumprod(s, i, j)
    prod = ones(length(s[1]))
    for k = i+1 : j-1
        prod .*= s[k]
    end
    return prod
end



""" terminal time commitor. 1 on inds and 0 on the rest for the terminal time """
function termcom(inds::Vector, n)
    c = zeros(n)
    c[inds] .= 1
    c
end

function commitor(T, terminal)
    TT = T
    nm = size(TT, 1)
    n = length(terminal)
    bnd = nm - n

    TT = TT - I
    for i in bnd+1 : nm
        TT[i,:] .= 0
        TT[i,i] = 1
    end
    x = zeros(nm)
    x[bnd+1:end] = terminal
    TT\x
end