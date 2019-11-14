using Arpack, LinearAlgebra
using SparseArrays

function galerkin(process, ts)
    n = length(process)
    m = length(ts) - 1

    dt = vcat(diff(ts))

    Q = [generatormatrix(process, t) for t in ts[1:m]]
    qout = [-diag(Q[t]) for t in 1:m]
    qt = [(Q[t] - Diagonal(Q[t])) ./ qout[t] for t in 1:m]
    s = [exp.(-dt[t]*qout[t]) for t in 1:m]


    T = spzeros(n*m, n*m)
    timeslice(i,j) = view(T, (i-1)*n + 1 : (i-1)*n + n, (j-1)*n + 1 : (j-1)*n + n)

    for ti in 1:m
        for tj = ti:m
            if ti==tj
                timeslice(ti,tj) .=
                    qt[tj] .* (1 ./ qout[ti] .* (s[ti] + dt[ti] * qout[ti] .- 1)) / dt[ti]
            elseif ti<tj
                timeslice(ti,tj) .=
                    qt[tj] .* (1 ./ qout[ti] .* (1 .-s[ti]) .* (1 .-s[tj]) .* dotprod(s, ti, tj)) / dt[ti]
            end
        end
    end

    if ts[end] == Inf  # account for absorbing boundary
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

function embeddedmatrix(q)
    qt = q - Diagonal(q)
    qout = sum(qt, dims=2)
    qt = qt ./ qout
    return qt, qout
end


# workaroud solution:
# iteratively select upper blocks
function sortschur(G)
    S = G

    N = length(S.values)
    sel = fill(true, N)
    for i = N:-1:1
        sel[1:i] .= true
        sel[i+1:N] .= false
        ind = argmin(real.(S.values[1:i]))

        sel[ind] = false
        if (ind<N) && isapprox(real(S.values[ind]), real(S.values[ind+1]))
            sel[ind+1] = false
        end
        S = ordschur!(S, sel)
    end
    S
end

function commitor(T, terminal)
    TT = T |> collect
    nm = size(TT, 1)
    n = length(terminal)
    bnd = nm - n

    TT[1:bnd, 1:bnd] = TT[1:bnd, 1:bnd] - I
    TT[bnd+1:end, :] .= 0
    TT[bnd+1:end, bnd+1:end] = I(n)
    x = zeros(nm)
    x[bnd+1:end] = terminal
    TT\x
end