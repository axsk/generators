import Arpack: eigs
using LinearAlgebra

function pccap(T::Matrix, n::Integer; opt = false)
    israte = isratematrix(T)
    pi     = stationarydensity(T, israte)
    X, λ   = schurvectors(T, pi, n, israte)
    chi    = makeprobabilistic(X, false)
end

function isratematrix(T::Matrix)
    s = sum(T[1,:])
    isapprox(s, 0) && return true
    isapprox(s, 1) && return false
    error("given matrix is neither a rate nor a probability matrix")
end

function stationarydensity(T, israte=isratematrix(T))
    which = israte ? :SM : :LM
    pi = eigs(T', nev=1, which=which)[2]
    @assert isreal(pi)
    pi = abs.(pi) |> vec
    pi = pi / sum(pi)
end

function makeprobabilistic(X::Matrix, opt::Bool)
    n = size(X, 2)
    A = innersimplexalgorithm(X)
    if n > 2 && opt
        A = opt(A, X, A->I3(A))
    end
    chi = X*A
end

function crispassignments(chi)
     assignments = mapslices(argmax, chi, dims=2) |> vec
end

function schurvectors(T, pi, n, israte)
    Tw = Diagonal(sqrt.(pi))*T*Diagonal(1 ./ sqrt.(pi)) # rescale to keep markov property
    Sw = schur!(Tw)                       # returns orthonormal vecs by def
    Xw, λ = selclusters!(Sw, n, israte)
    X  = Diagonal(1 ./sqrt.(pi)) * Xw              # scale back
    X  = X[1,1]>0 ? X : -X
    X, λ
end

# select the schurvectors corresponding to the n abs-largest eigenvalues
# if reverse==true select highest abs value, otherwise select lowest (for rate matrices)
function selclusters!(S, n, ratematrix)
    ind = sortperm(abs.(S.values), rev=!ratematrix) # get indices for dominant eigenvalues
    select = zeros(Bool, size(ind))           # create selection vector
    select[ind[1:n]] .= true
    S = ordschur!(S, select)                  # reorder selected vectors to the left
    if !isapprox(S.T[n+1, n], 0)                       # check if we are cutting along a schur block
        @error("conjugated eigenvector missing")
        display(S.T)
    end
    S.vectors[:,1:n], S.values[1:n]       # select first n vectors
end

# compute initial guess based on indexmap
innersimplexalgorithm(X) = feasiblize!(inv(X[indexmap(X), :]), X)

function indexmap(X)
    # get indices of rows of X to span the largest simplex
    rnorm(x) = sqrt.(sum(abs2.(x), dims=2)) |> vec
    ind=zeros(Int, size(X,2))
    for j in 1:length(ind)
        rownorm=rnorm(X)
        # store largest row index
        ind[j]=argmax(rownorm)
        if j == 1
            # translate to origin
            X = X .- X[ind[1],:]'
        else
            # remove subspace
            X=X/rownorm[ind[j]]
            vt=X[ind[j],:]'
            X=X-X*vt'*vt
        end
    end
    return ind
end

function feasiblize!(A,X)
    A[:,1] = -sum(A[:,2:end], dims=2)
    A[1,:] = -minimum(X[:,2:end] * A[2:end,:], dims=1)
    A / sum(A[1,:])
end

# crispness criterion, cf. Roeblitz (2013)
# only applies if X is normalized (eq. 8)
# TODO: check does this apply?
function roeblitzcrit(A)
    n = size(A,1)
    trace = 0
    for i = 1:n, j=1:n
        trace += A[i,j]^2 / A[1,j]
    end
    return trace
end

function opt(A0, X)
    A = copy(A0)
    Av = view(A, 2:size(A,1), 2:size(A,2)) # view on the variable part

    function obj(a)
        Av[:] = a
        -roeblitzcrit(feasiblize!(A, X))
    end

    result = optimize(obj, Av[:], NelderMead())
    Av[:] = result.minimizer
    return feasiblize!(A, X)
end

function randomstochasticmatrix(n, reversible=true)
    P = rand(n,n)
    if reversible
        P = (P + P') / 2
    end
    P ./= sum(P, dims=2)
end

function isreversible(P, pi=stationarydensity(P))
    db = [pi[i] * P[i,j] - pi[j] * P[j,i] for i =1:5, j=1:5]
    all(isapprox.(db, 0, atol=1e-8))
end

function test()
    χ = pccap(randomstochasticmatrix(8), 2)
    a = crispassignments(χ)
    χ, a
end

test()
