using SparseArrays

""" Given a galerkin matrix of size (nm)x(nm), compute the elementary backward commitors for the (n) space-points over the (m) time spans
Returns a matrix whose `n` columns are the (nm) space-time entries of the `n` commitors """
function commitors(g, n)
    C = g - I
    nm = size(C, 1)
    bnd = nm - n

    # fix non-transitions by assigning them to the target set
    s = sum(C, dims=2)
    for i in 1 : nm
        C[i, nm - n + (i-1)%n + 1] = -s[i]
    end

    # set boundary condition
    for i in bnd+1 : nm
        C[i,:] .= 0
        C[i,i] = 1
    end

    term = spzeros(nm,n)
    term[nm-n+1:end, :] = sparse(I, n,n)

    cs = C\collect(term)
end