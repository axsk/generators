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

length(b::BickleyJet) = b.nx * b.ny


function generator(t, conf::BickleyJet)
    @unpack_BickleyJet conf

    psi(t,x,y) = -U0 * L * tanh(y/L) + U0 * L * sech(y/L)^2 *
        (A2 * cos(k2*(x-c2*t)) + A3 * cos(k3*(x-c3*t)))

    xs = range(xmin, stop=xmax, length=nx+1)
    ys = range(ymax, stop=ymin, length=ny+1)

    dx = step(xs)
    dy = step(ys)

    # optimized functions for the flow using independence on other variables
    @memoize function vflow(i,j,t)
        xx = (xs[j]+xs[j+1]) / 2
        flow = derivative(x->psi(t,x,ys[i+1]), xx) * abs(dx)
        #@show "v",t,i,j,flow
        pos = max( flow, 0)
        neg = max(-flow, 0)
        pos, neg
    end

    @memoize function hflow(i,j,t)
        yy = (ys[i]+ys[i+1]) / 2
        flow = derivative(y->-psi(t,xs[j+1],y), yy) * abs(dy)
        #@show "h",t,i,j,flow
        pos = max( flow, 0)
        neg = max(-flow, 0)
        pos, neg
    end

    G = spzeros(ny*nx, ny*nx)

    for i = 1:ny-1
        for j = 1:nx-1
            # from anchor i,j to bottom and right
            cell = i + (j-1)*ny
            pos, neg = vflow(i,j,t)
            G[cell, cell+1]  = neg
            G[cell+1, cell]  = pos

            pos, neg = hflow(i,j,t)
            G[cell, cell+ny]  = pos
            G[cell+ny, cell]  = neg
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
