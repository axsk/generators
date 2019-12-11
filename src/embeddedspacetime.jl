using Arpack, LinearAlgebra
using SparseArrays

export galerkin

function galerkin(Qs, dt)
    n = size(Qs[1], 1)
    m = length(Qs)

    @assert length(dt) == m

    qout = [-diag(Qs[t]) for t in 1:m]
    qt = [(Qs[t] - Diagonal(Qs[t])) ./ qout[t] for t in 1:m]
    s = [exp.(-dt[t]*qout[t]) for t in 1:m]

    T = spzeros(n*m, n*m)
    timeslice(i,j) = view(T, (i-1)*n + 1 : (i-1)*n + n, (j-1)*n + 1 : (j-1)*n + n)

    for ti in 1:m
        for tj in ti:m
            if ti==tj
                timeslice(ti,tj) .=
                    qt[tj] .* (1 ./ qout[ti] .* (s[ti] + dt[ti] * qout[ti] .- 1)) / dt[ti]
            elseif ti<tj
                timeslice(ti,tj) .=
                    qt[tj] .* (1 ./ qout[ti] .* (1 .-s[ti]) .* (1 .-s[tj]) .* dotprod(s, ti, tj)) / dt[ti]
            end
        end
    end

    if dt[end] == Inf  # account for absorbing boundary
        timeslice(m,m) .= qt[m]
    end
    T
end

function dotprod(s, i, j)
    prod = ones(length(s[1]))
    for k = i+1 : j-1
        prod .*= s[k]
    end
    return prod
end



""" terminal time commitor. 1 on inds and 0 on the rest for the terminal time """
function termcom(inds, n)
    c = zeros(n)
    c[inds] .= 1
    c
end

function commitor(T, terminal)
    TT = T |> collect
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