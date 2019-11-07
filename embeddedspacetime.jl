using Arpack, LinearAlgebra
using SparseArrays
using Plots
# actually it should be embedded augmented / space-time markov chain right?
# embedded space-time jump chain (estjc)

function augmentedembeddedmatrix(process, ts)
    ntime = length(ts)
    n  = length(process)

    dt = vcat(diff(ts), Inf) # boundary open to the right
    T = spzeros(n*ntime, n*ntime)

    function timeslice(i,j)
        view(T, (i-1)*n + 1 : (i-1)*n + n, (j-1)*n + 1 : (j-1)*n + n)
    end

    Q = [generatormatrix(process, t) for t in ts]
    H = [Diagonal(exp.(dt[i] * diag(Q[i]))) for i in 1:ntime] # prob not to jump (hold) in timeframe 
    L = [I - H for H in H] # prob to jump (leave) inside timeframe
    L[end] = Diagonal(ones(n))

    R = [I - Diagonal(Q)\Q for Q in Q] # jump matrix

    for i = 1:ntime
        for j = i:ntime
            if i == j
                timeslice(i,j) .= L[i] * R[i]
            else
                timeslice(i,j) .= L[j] * prod(H[k] for k=i:j-1) * R[j]
            end
        end
    end
    T
end

function galerkin(process, ts)
    n = length(process)
    m = length(ts) - 1

    dt = vcat(diff(ts))

    Q = [generatormatrix(process, t) for t in ts]
    qout = [-diag(Q[t]) for t in 1:m]
    qt = [(Q[t] - Diagonal(Q[t])) ./ qout[t] for t in 1:m]
    s = [exp.(-dt[t]*qout[t]) for t in 1:m]


    T = spzeros(n*m, n*m)
    timeslice(i,j) = view(T, (i-1)*n + 1 : (i-1)*n + n, (j-1)*n + 1 : (j-1)*n + n)

    for ti in 1:m
        for tj = ti:m
            if ti==tj
                timeslice(ti,tj) .=
                    qt[tj] .* (1 ./ qout[ti] .* (s[ti] + dt[ti] * qout[ti] .- 1))
            elseif ti<tj
                timeslice(ti,tj) .=
                    qt[tj] .* (1 ./ qout[ti] .* (1 .-s[ti]) .* (1 .-s[tj]) .* dotprod(s, ti, tj))
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
