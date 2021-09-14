using Arpack, LinearAlgebra
using SparseArrays

export galerkin

galerkin(Qs::Vector{<:Matrix}, dt) = galerkin(map(sparse, Qs), dt)
function galerkin(Qs::Vector{<:SparseMatrixCSC}, ts)
    dt = diff(ts)
    n = size(Qs[1], 1)
    m = length(dt)


    I = Int[]
    J = Int[]
    V = Float64[]

    qout = [collect(-diag(Qs[t])) for t in 1:m]
    s = [exp.(-dt[t]*qout[t]) for t in 1:m] # note that we still need the possible 0 rates here (c.f. below). in case of dt == Inf this can lead to 0*Inf = NaN
    qt = [dropzeros(replace((Qs[t] - Diagonal(Qs[t])) ./ qout[t], NaN => 0)) for t in 1:m]

    for ti in 1:m
        for tj in ti:m
            if ti==tj
                fact = (s[ti] + dt[ti] * qout[ti] .- 1) ./ qout[ti]
                # if qout == 0, fact should calculate (1-1)/0 = NaN, but actually we have no transition, so set it to 0
                replace!(fact, NaN => 0)
            elseif ti<tj
                t = (1 .- s[ti]) ./ qout[ti]
                # if qout == 0, t should calculate NaN, but actually the integral solves to 1
                replace!(t, NaN => 1)
                fact = t .* (1 .- s[tj]) .* mycumprod(s, ti, tj)
            end
            res  = qt[tj] .* (fact ./ dt[ti])
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