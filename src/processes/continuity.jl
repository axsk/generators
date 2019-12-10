# processes governing the flow of incompressible fluids
# here we have dp/dt + grad p * u + p * grad u = 0
# with density p and flow field u
# for incompressible flows grad u = 0

Base.range(a,b) = range(a,b, step=1)
Base.range(a,b,n::Int) = range(a, b, length=n)

@with_kw struct Flow <: Process
    u = (t, x, y)->[1,0]
    xs = range(0,nx)
    ys = range(0,ny)
end

generatormatrix(f::Flow, t) = flowgenerator(f.xs, f.ys, (x,y) -> f.u(x,y,t))

function flowgenerator(xs, ys, u)
    nx = length(xs)
    ny = length(ys)

    G = spzeros(nx*ny, nx*ny)
    ind = LinearIndices((1:nx, 1:ny))
    for (i,x) in enumerate(xs)
        for (j,y) in enumerate(ys)
            ux, uy = u(x,y)
            ii = i + Int(sign(ux)) 
            if 1 <= ii <= nx 
                G[ind[i,j], ind[ii,j]] = abs(ux)
            end
            jj = j + Int(sign(uy))
            if 1 <= jj <= ny 
                G[ind[i,j], ind[i,jj]] = abs(uy)
            end
        end
    end
    for i = 1:size(G,1)
        G[i,i] = - sum(G[i,:])
    end
    G
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

function myquiver(f::Flow, t)
    myquiver(f.xs, f.ys, (x,y)->f.u(t,x,y))
end
