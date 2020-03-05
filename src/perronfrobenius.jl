# Reconstruction of the Perron Frobenius operator
# by means of the activity map

using LinearAlgebra

const neumann_start = 1
const NEUMANN_N = 100

# neumann series
function jumpactivity(g)
    power = one(g)
    j = zero(g)
    for i in 1:NEUMANN_N
        power = power * g
        j += power
    end
    j
end

function project_jumpactivity(E, f0)
    nx = length(f0)
    nt = round(Int, size(E, 1)/nx)
    activity = sum(E[1:nx, :] .* f0, dims=1)
    reshape(activity, nx, nt)
end

function perronfrobenius(g, qs, ts)
    nx = size(qs[1],1)
    nt = length(ts) - 1
    qout = [collect(-diag(q)) for q in qs]

    # survival probability from one galerkin cell to T=ts[end]
    S = zeros(nx, nt)
    for i=1:nt
        # int_t[i]^t[i+1] exp(-qT + qx) dx
        S[:,i] .= 1 ./ qout[i] .*(exp.(-qout[i]*(ts[end] - ts[i+1])) - exp.(-qout[i]*(ts[end] - ts[i]))) / (ts[i+1]-ts[i])
    end

    # survival probaility from ts[1] to ts[end]
    surviveall = Diagonal(exp.(-qout[1]*(ts[end] - ts[1])))

    E = jumpactivity(g)

    Pf = sum(reshape(E, nx, nt, nx, nt)[:,1,:,:] .* reshape(S, 1, nx, nt), dims=3) |> x->reshape(x, nx, nx)
    Pf = Pf + surviveall
    Pf' # since we multiply from right
end