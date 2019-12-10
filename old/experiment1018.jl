include("main.jl")

tswitch = [1,2,2.5,3,3.2, 3.3]
cellres = 2
rate = 0.5

p = BarrierSwitch(tswitch, cellres, rate)

ts = 0:.1:4

T = augmentedembeddedmatrix(p, ts)

using Makie
Makie.heatmap(T::SparseMatrixCSC) = hm(Matrix(T[end:-1:1, :]') |> Matrix)

hm = Makie.heatmap

hm(log.(sum(T^n for n=0:100)))