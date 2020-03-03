using SparseArrays

""" Given a galerkin matrix of size (nm)x(nm), compute the elementary backward commitors for the (n) space-points over the (m) time spans
Returns a matrix whose `n` columns are the (nm+n) space-time entries of the `n` commitors """
function commitors(g, n)
    nm = size(g, 1)
    A = spzeros(nm+n, nm+n)
    A[1:nm, 1:nm] = g-I # gx = x <=> (g-I)x = 0

    # assing non-transitions to the target set
    s = sum(A, dims=2)
    for i in 1 : nm
        A[i, nm + (i-1)%n + 1] = -s[i]
    end

    b = spzeros(nm+n,n)
    # set boundary condition
    for i in 1:n
        A[nm+i,nm+i] = 1
        b[nm+i,i]    = 1
    end

    b[nm+1:end, :] .= sparse(I, n,n)

    # replace the solve by own solver
    #cs = A\collect(b) # we need dense b, since Julia does not yet support sparse RHS

    m = Int(nm/n)
    mysolve(A, n, m+1)
end

# solve Ax=b for x where b has the terminal structure (Id on the bottom and 0 else)
function mysolve(A, n,m)
    # split A into timeslices
    As = [A[i*n .+ (1:n), j*n .+ (1:n)] for i=0:m-1, j=0:m-1]
    xs = Array{Matrix{Float64}}(undef, m)
    xs[end] = Matrix(I, n, n) # boundary condition
    for i=m-1:-1:1
        b = zeros(n,n)
        for j=i+1:m
            b -= As[i,j] * xs[j]
        end
        xs[i] = As[i,i] \ b
    end
    vcat(xs...)
end