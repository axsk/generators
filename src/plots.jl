#using Makie
import Makie
import Makie: wireframe, wireframe!, lift

function visualize(r)
    s = Makie.Scene()
    xs = r.xs |> collect
    ts = r.ts[1:end-1]
    V = [r.V(x,t) for x in xs, t in ts]
    Makie.wireframe!(s, xs, ts, V)
    C = reshape(r.c, length(xs), length(ts))
    Makie.wireframe!(s, xs, ts, C, color=:blue)
    s[Axis][:names][:axisnames] = ("x", "t", "f")
    s
end

function interact()
    s = Makie.Scene()
    sl = slider(1:10, raw = true, camera = campixel!, start = 8)
    C = lift(sl[end][:value]) do i
        com = zeros(10) 
        com[i] = 1
        r = experiment(PotentialShiftingBarrier(), grid, ts, com)
        reshape(r.c, length(grid), length(ts)-1)
    end
    hbox(Makie.wireframe(grid|>collect, ts[1:end-1], C), sl)
end

import Makie.AbstractPlotting.plot!

wireframe(a::AbstractVector, b::AbstractVector, c::Function) = wireframe(a, b, [c(a,b) for a in a, b in b])
wireframe(a::AbstractVector, b::AbstractVector, c::AbstractVector) = wireframe(a, b, reshape(c, length(a), length(b)))
wireframe(c::Matrix) = wireframe(1:size(c,1), 1:size(c,2), c)

function plot!(plot::Makie.Wireframe{ <: Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}})
    c = lift((a,b,c)->reshape(c, length(a), length(b)), plot[1], plot[2], plot[3])
    wireframe!(plot, plot[1], plot[2], c)
end

function plot!(plot::Makie.Wireframe{ <: Tuple{Array{Float64,1},Array{Float64,1},Function}})
    c = lift((a,b,c)->[c(a,b) for a in a, b in b], plot[1], plot[2], plot[3])
    wireframe!(plot, plot[1], plot[2], c)
end

function plot!(plot::Makie.Wireframe{ <: Tuple{Array{Float64,2}}})
    a = lift(x->1:size(x,1), plot[1])
    b = lift(x->1:size(x,2), plot[1])
    wireframe!(plot, a, b, plot[1])
end

function plot!(plot::Makie.Wireframe{ <: Tuple{Array{Float64,3}}})
    c = lift(plot[1]) do x
        for i=1:size(x, 3)
            wireframe!(plot, x[:,:,i])
        end
    end
end

function myquiver(xs, ys, u)
    tx, ty, tux, tuy = [], [], [], []
    for x in xs
        for y in ys
            push!(tx, x)
            push!(ty, y)
            ux, uy = u(x,y)
            push!(tux, ux)
            push!(tuy, uy)
        end
    end
    quiver(tx,ty,quiver=(tux,tuy))
end