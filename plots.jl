using Makie

function spacetimevector(x, n, m)
    @assert length(x) == n * m
    sc = Scene()

    for t = 1:m
        xx = x[(t-1)*n.+(1:n)]
        lines!(sc, 1:n, xx, ones(n) * t)    
    end
    sc
end