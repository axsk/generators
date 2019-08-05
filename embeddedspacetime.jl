
using Arpack, LinearAlgebra
using Plots
# actually it should be embedded augmented / space-time markov chain right?

function augmentedembeddedmatrix(process, ts)
    ntime = length(ts)
    n  = spacesize(process)

    dt = vcat(diff(ts), Inf)
    T = zeros(n*ntime, n*ntime)

    function timeslice(i,j)
        view(T, (i-1)*n + 1 : (i-1)*n + n, (j-1)*n + 1 : (j-1)*n + n)
    end

    Q = [generatormatrix(process, t) for t in ts]
    H = [Diagonal(exp.(dt[i] * diag(Q[i]))) for i in 1:ntime]
    L = [I - H for H in H]
    L[end] = Diagonal(ones(n))

    R = [I - Diagonal(Q)\Q for Q in Q]

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

using Parameters

@with_kw struct DoubleWell
    beta = 1
    ngrid = 20
end

spacesize(d::DoubleWell) = d.ngrid

function generatormatrix(d::DoubleWell, t)
    @unpack beta, ngrid = d
    n = ngrid

    V(x) = (x^2-1)^2
    grid = range(-2,2, length=ngrid)
    VV = map(x->exp(-beta * V(x)), grid)
    G = zeros(n,n)
    for i = 1:n-1
        G[i, i+1] = sqrt(VV[i+1]   / VV[i])
        G[i+1, i] = sqrt(VV[i] / VV[i+1])
    end
    G[1,end] = sqrt(VV[end] / VV[1])
    G[end,1] = sqrt(VV[1]   / VV[end])
    for i = 1:n
        G[i,i] = -sum(G[i,:])
    end
    G
end

using Arpack

function myeigs(A; nev=5, which = :LR, kwargs...)
    e, v = eigs(A; nev=nev, which=which, kwargs...)
    #e = real(e)
    #v = real(v)

    for i = 1:size(v,2)
        v[:,i] = v[:,i] / v[findmax(abs.(v[:,i]))[2], i]
    end
    e,v
end

function ploteigs(e,v)
    e = real(e)
    v = real(v)
    plot(v, line_z=e', color=:darktest)
end

function sortschur(G, n)
    S = schur(collect(G))
    ind = sortperm(real(S.values), rev=true) # get indices for dominant eigenvalues
    select = zeros(Bool, size(ind))           # create selection vector
    select[ind[1:n]] .= true
    S = ordschur!(S, select)
     S.values[1:n], S.vectors[:,1:n]
end

#function test()
process = DoubleWell(beta = 3)
timegrid = range(0, step=1, length=3)
G = augmentedembeddedmatrix(process, timegrid)
heatmap(G)
#e, v = myeigs(G, nev=9)
#e2 = e[2][:,2] |> real
#ploteigs(e,v)
e,v = sortschur(G, 9)
@show e
plot(v)

#end

#test()I
