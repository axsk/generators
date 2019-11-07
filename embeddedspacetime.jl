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
