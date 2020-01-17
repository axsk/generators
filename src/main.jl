using Parameters
#include("embeddedspacetime.jl")
#include("plots.jl")

gauss(x, mu, sigma) = 1 / (sigma * sqrt(2*pi)) * exp( - 1/2 * ((x-mu)/sigma)^2)
PotentialShiftingBarrier(sigma) = (x,t) -> gauss(x, 0.2 + t/2, sigma)

PotentialDoubleWell() = (x,t) -> (x^2-1)^2 
#ShiftingWells() = (x,t) -> (x^2-1)^2+x*t 

AppearingBarrier() = (x,t) -> t*gauss(x, .5, 0.1)
VanishingBarrier(;mu=.5, sigma=.1) = (x,t) -> (1-t) * gauss(x, mu, sigma)

@with_kw struct Experiment
    xs = range(-1.5, 1.5, length=20) |> collect
    ts = 0:.1:2
    setA = [18]
    V = PotentialShiftingBarrier(1)
    beta = 1
    flux = 1

    Qs =[sqra(x->V(x,t), xs, beta, flux) for t in ts] 
    dt = vcat(diff(ts), Inf)
    g = galerkin(Qs, dt)
    c = commitor(g, termcom(setA, length(xs)))
    cs = hcat((commitor(g, termcom([i], length(xs))) for i in 1:length(xs))...)
end

n(e::Experiment) = length(e.xs)
m(e::Experiment) = length(e.ts)

cplot(e::Experiment) = wireframe!(e.xs, e.ts, e.c)

import ColorSchemes
ncolors(n) = get(ColorSchemes.plasma, 1:n, :extrema)

function viscoms()
    s = Makie.Scene()
    cs = ncolors(n(e))
    for i = 1:n(e)
        wireframe!(s, e.xs|>collect, e.ts|>collect, reshape(e.cs[i], n(e), m(e)), color=cs[i])
    end
    s
end

using JuMP

function optim(cs, n, obj=:offdiag)

    l = size(cs, 2)
    m=Model()
    @variable(m, x[1:l,1:n])
    #@objective(m, Max, sum(diag((coms * x)' * coms * x)))
    #@objective(m, Max, sum((cs*x).^2))
    @expression(m, z, cs * x)
    if obj == :offdiag
        @objective(m, Min, sum(z[:,i]' * z[:,j] for i=1:n, j=i+1:n))
    elseif obj == :diag
        @objective(m, Max, sum(z.^2))
        @constraint(m, mass, sum(x, dims=1) .>= 0.1)
    elseif obj == :det
        @objective(m, Max, det(z'*z))
    end
    @constraint(m, c1, x .>=0)
    @constraint(m, c2, x .<=1)
    @constraint(m, partition, sum(x, dims=2) .== 1)
    m
end

function nmobj(c, A)
    x = c*A
    return sum(maximum(c*A[:,j]) for j=1:size(x,2))
end

function feas(A)
    A = A .- min.(0, minimum(A, dims=2)) # shift to positive
    A = A ./ sum(A, dims=2) # rescale to rowsum 1
end

using Optim
function nm(C, m)
    n = size(C, 1)
    A = rand(n,m)
    lower = zeros(n)
    upper = ones(n)
    sol = optimize(x->-nmobj(C, feas(reshape(x, n, m))), vec(A), NelderMead(), Optim.Options(iterations=100000))
    feas(reshape(sol.minimizer, n, m)) , sol
end

function totalcorrelation(X)
    let n = size(X, 1)
      X= X ./ sum(X)
      M = X ./ sum(X, dims=1)
      c = 0
      for x1 = 1:n
          for x2 = 1:n
              p = X[x1, 1] * X[x2, 2]
              c += p * log(p / (M[x1, 1] * M[x2,2]))
          end
      end
  end
    c
end

function maxindices(c, m)
  
  " "
  n = size(c, 1)
end

function vertexindices(X, m)
  "get indices of rows of X to span the largest simplex"
  ind = zeros(Int, m)
  for j in 1:m
      rownorm=norm.(eachrow(X))
      # store largest row index
      ind[j]=argmax(rownorm)
      if j == 1
          # translate to origin
          X = X .- X[ind[1],:]'
      else
          # remove subspace
          X=X/rownorm[ind[j]]
          vt=X[ind[j],:]'
          X=X-(X*vt')*vt
      end
  end
  return ind
end


function optimize_maxassignment(c, m)
  "return the A matrix maximizing the objective o(A) = sum_i=1^m max(chi(:, i))
   this means that the commitors are distributed such that the sum of the maximal membership of each cluster is maximized

  ind is a heuristic list[m] of the indicies at which the succesive clusters are maximal
  each commitor (column of c) gets assigned to one of the clusters
  assigning commitor c[:, i] to cluster m, the growth of the objective is c[ind[m], m].
  so we want to assign it to the cluster m such that that term is maximized for each c_i"
  
  n = size(c, 1)
  A = zeros(n, m)

  ind = vertexindices(c, m)

  for i in 1:n
    m = argmax(c[ind, i])
    A[i, m] = 1
  end

  return A
end
