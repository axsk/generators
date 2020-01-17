using Statistics
using SparseArrays
using Plots
using ForwardDiff: derivative
using Parameters
using Memoize

@with_kw struct BickleyJet
    xmin = - pi * 6.371
    xmax = pi * 6.371
    ymin = -3.
    ymax = 3.

    nx = 120
    ny = 36

    U0 = 5.4138
    L = 1.77
    c2 = 0.205 * U0
    c3 = 0.461 * U0
    A2 = 0.1
    A3 = 0.3
    k2 = 4 / 6.371
    k3 = 6 / 6.371
end

xs(b::BickleyJet) = range(b.xmin, stop=b.xmax, length=b.nx) 
ys(b::BickleyJet) = range(b.ymin, stop=b.ymax, length=b.ny)

length(b::BickleyJet) = b.nx * b.ny

function bickleyflow(b::BickleyJet, t, x, y)
    @unpack_BickleyJet b

    psi(t,x,y) = -U0 * L * tanh(y/L) + U0 * L * sech(y/L)^2 *
        (A2 * cos(k2*(x-c2*t)) + A3 * cos(k3*(x-c3*t)))
    
    # TODO: scale with cell size
    dx = derivative(y->-psi(t,x,y), y)
    dy = derivative(x-> psi(t,x,y), x)

    (dx, dy)
end



function generator(t, conf::BickleyJet)
    @unpack_BickleyJet conf

    psi(t,x,y) = -U0 * L * tanh(y/L) + U0 * L * sech(y/L)^2 *
        (A2 * cos(k2*(x-c2*t)) + A3 * cos(k3*(x-c3*t)))

    xs = range(xmin, stop=xmax, length=nx+1)
    ys = range(ymax, stop=ymin, length=ny+1)

    dx = step(xs)
    dy = step(ys)

    # indices xs[i], ys[j] denote the top left corner of the box
    # vx(i,j,t) is the horizontal velocity at the center of the edge between box (i,j) and (i,j+1), i.e. at x[j+1]

    # optimized functions for the flow using independence on other variables
    vy(x,y) = derivative(x-> psi(t,x,y), x)
    vx(x,y) = derivative(y->-psi(t,x,y), y)

    ind = LinearIndices((1:nx, 1:ny))
    ind = vcat(ind, ind[1,:]') # add indices for periodic x direction
    G = spzeros(ny*nx, ny*nx)

    for i = 1:nx
        for j = 1:ny
            # from anchor i,j to right
            from = ind[i,j]
            to   = ind[i+1, j]
            v = vx(xs[i+1], (ys[j]+ys[j+1])/2)
            if v < 0
                from, to = to, from
            end
            G[from, to] = abs(v) * abs(dx)

            # and to the bottom
            j == ny && continue # reached the bottom edge
            from = ind[i,j]
            to   = ind[i,j+1]
            v = vy((xs[i]+xs[i+1])/2, ys[j+1])
            if v < 0 
                from, to = to,from
            end
            G[from, to] = abs(v) * abs(dy)
        end
    end

    # fix diagonal
    for i = 1:nx*ny
        G[i,i] = -sum(G[i,:])
    end
    G
end

using KrylovKit

function stepflow(conf::BickleyJet, x0, t0, dt)
    @unpack nx, ny = conf
    g = generator(t0, conf)
    xt, _ = exponentiate(g', dt, vec(x0))
    xt = reshape(xt, ny, nx)
end

x0 = zeros(36,120); x0[1:6:end,1:2:end].=1; x0[1:2:end, 1:20:end] .= 1

function testflow(conf::BickleyJet = BickleyJet();
    x0 = sprand(conf.ny, conf.nx, 0.1),
    t0 = 0,
    dt = 0.1,
    n = 10)
    x = [x0]
    for i=1:n
        xt = stepflow(conf, x[end], t0, dt)
        heatmap(xt, title="$i") |> display
        push!(x, xt)
        t0 = t0 + dt
    end

    x
end

function flowheatmap(generator::AbstractMatrix, nx=120, ny=36)
    heatmap(reshape(diag(G), ny, nx))
end
