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

function PotentialExperiment(;
    xmin = -1,
    xmax = 1,
    n = 20,
    xs = range(xmin, xmax, length=n) |> collect,
    ts = 0:.1:2,
    V = AppearingBarrier(),
    beta = 1,
    flux = 1,
    Qs = [sqra(x->V(x,t), xs, beta, flux) for t in ts],
    dt = vcat(diff(ts), Inf),
    g = galerkin(Qs, dt),
    c = hcat((commitor(g, termcom([i], length(xs))) for i in 1:length(xs))...),
    cc = c * optimize_maxassignment(c[1:length(xs),:], 2)[:,2] |> x->reshape(x, length(xs), length(ts)))

    d = Base.@locals
    e = NamedTuple{Tuple(keys(d))}(values(d))
    plot_pot(e) |> display
    e
end

function plot_pot(e)
    plot(e.xs, e.cc, title = "n=$(e.n) b=$(e.beta) f=$(e.flux) V=l$(methods(e.V).ms[1].line)", labels=e.ts')
end