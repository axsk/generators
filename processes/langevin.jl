using Parameters

@with_kw struct OverdampedLangevin <: Process
    V::Function
    phi
    grid::Vector
    periodic::Bool = false
end

length(p::OverdampedLangevin) = length(p.grid)


# todo: add temperature/noise and flux-constant!
function generatormatrix(p::OverdampedLangevin, t)
    @unpack_OverdampedLangevin p

    n = length(grid)
    VV = map(x->exp(-V(x,t)), grid)

    G = zeros(n,n)

    # square root approximation
    for i = 1:n-1
        G[i, i+1] = sqrt(VV[i+1]   / VV[i])
        G[i+1, i] = sqrt(VV[i] / VV[i+1])
    end

    if periodic
        G[1,end] = sqrt(VV[end] / VV[1])
        G[end,1] = sqrt(VV[1]   / VV[end])
    end

    # diagonal
    for i = 1:n
        G[i,i] = -sum(G[i,:])
    end
    G * phi
end

# autonomous double well
# note that here beta is used to scale the potential
DoubleWell(beta, phi, ngrid) = OverdampedLangevin(
    (x,t) -> beta * (x^2-1)^2,
    phi,
    range(-2,2, length=ngrid),
    false)

ShiftingWells() = OverdampedLangevin(
    (x,t) -> (x^2-1)^2+x*t,
    1,
    range(-2,2, length=5),
    false)

function RisingWell(; beta=1, phi=1, speed=1, n=6, lim=2)
    OverdampedLangevin(
        (x,t) -> beta * (x^2-1)^2 + speed * t * x,
        phi,
        range(-lim,lim, length=n),
        false)
end