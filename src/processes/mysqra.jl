using LinearAlgebra

# convenience functions
linspace(a, b, n) = range(a, stop=b, length=n)
unitdoublewell(x) = (x^2-1)^2 # double well with minima at +-1
makeQ(Q) = Q - Diagonal(sum(Q, dims=2)|>vec)
statDens(Q) = eigen(Q').vectors[:,end]

# setup
n  = 100
xs = linspace(-2,2,n)
U = [unitdoublewell(x) for x in xs]
beta = 1


pi = exp.(-beta * U)
Q_sqra = [ abs(i-j)!=1 ? 0 : sqrt(pi[j]/pi[i]) for i=1:n, j=1:n] |> makeQ
Q_adv  = [ abs(i-j)!=1 ? 0 : sqrt(pi[j]/pi[i]) * (U[j]-U[i])/step(xs) for i=1:n, j=1:n] .|> abs |> makeQ
Q_diff = [ abs(i-j)!=1 ? 0 : sqrt(pi[j]/pi[i]) *  for i=1:n, j=1:n] |> makeQ

Q = Q_gal + Q_sqra

q = statDens(Q) ./ pi