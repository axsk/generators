rmnan(X) = map(x->isnan(x) ? 0 : x, X)

function fixgalerkin(G)
    G = rmnan(G)
    for i = 1:size(G,1)
        G[i,i] = 1 - sum(G[i, :])
    end
    G
end

function BickleyExperiment()
    b = BickleyJet(nx=20, ny=6)
    f = Flow(u = (t,x,y) -> bickleyflow(b, t, x, y), xs = xs(b), ys=ys(b))

    # either or
    q(t) = generator(t, b) # compute flow at boundaries
    #q(t) = generatormatrix(f, t) # compute flow at midpoints

    ts = [0,1,2]
    dt = [1,1,1]

    qs = [q(t) for t in ts]
    GG = galerkin(qs, dt)

    G = fixgalerkin(GG)
    Base.@locals
end


