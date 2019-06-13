using Statistics

function psi(t,x,y)
    U0 = 5.4138
    L = 1.77
    c2 = 0.205 * U0
    c3 = 0.461 * U0
    A2 = 0.1
    A3 = 0.3
    k2 = 4 / 6.371
    k3 = 6 / 6.371

    p = -U0 * L * tanh(y/L) + U0 * L * sech(y/L)^2 * (
        A2 * cos(k2*(x-c2*t)) + A3 * cos(k3*(x-c3*t)))
end

using ForwardDiff: derivative

const xmin = - pi * 6.371
const xmax = pi * 6.371
const ymin = -3.
const ymax = 3.

const nx = 120
const ny = 36

const xs = (range(xmin, stop=xmax, length=nx+1))
const ys = (range(ymax, stop=ymin, length=ny+1))

dx = step(xs)
dy = step(ys)

function generator(t)
    G = zeros(ny*nx, ny*nx)

    for i = 1:ny-1
        for j = 1:nx-1
            # from anchor i,j to bottom and right
            cell = i + (j-1)*ny
            @show pos, neg = vflow(i,t)
            G[cell, cell+1]  = neg
            G[cell+1, cell]  = pos

            @show pos, neg = hflow(j,t)
            G[cell, cell+ny]  = pos
            G[cell+ny, cell]  = neg
        end
    end
    # todo: fix diagonal
    for i = 1:nx*ny
        G[i,i] = -sum(G[i,:])
    end
    G
end

const nquad = 100

dpsidx(t,x,y) = derivative(x->psi(t,x,y), x)
dpsidy(t,x,y) = derivative(y->psi(t,x,y), y)

using Memoize

# optimized functions for the flow using independence on other variables
@memoize function vflow(i,t)
    flow = derivative(y->psi(t,0,y), ys[i+1])
    scale = abs(dy)
    pos = max(flow, 0) * scale
    neg = max(-flow, 0) * scale
    pos, neg
end

@memoize function hflow(j,t)
    flow = derivative(x->psi(t,x,0), xs[j+1])
    scale = abs(dx)
    pos = max(flow, 0) * scale
    neg = max(-flow, 0) * scale
    pos, neg
end

#=
# general functions (todo: integrate using mean)
function vflow(i,j,t)
    y = ys[i+1]
    (x0, x1) = xs[j], xs[j+1]
    flow = dy(t,0,y)
    scale = abs(x0 - x1)
    pos = max(flow, 0) * scale
    neg = max(-flow, 0) * scale
    pos, neg
end

function hflow(i,j,t)
    x = xs[j+1]
    (y0, y1) = ys[i], ys[i+1]
    flow = dx(t,x,0)
    scale = abs(y0 - y1)
    pos = max(flow, 0) * scale
    neg = max(-flow, 0) * scale
    pos, neg
end
=#

function transfer(rng)
    T = diagm(0=>ones(nx*ny))
    for t in rng
        T = T * exp(step(rng) * generator(t))
    end
    T
end
function transfer(x0, rng)
    x = vec(x0)'
    xhist = []
    for t in rng
        x = x * exp(step(rng) * generator(t))
        push!(xhist, reshape(x, ny, nx))
    end
    xhist
end

using SparseArrays
using Plots
x = sprand(ny, nx, 0.1)
g = generator(0)

#reshape(vec(x) * exp(g), ny, nx) |> collect |> heatmap
