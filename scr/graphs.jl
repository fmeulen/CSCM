using Laplacians
using PyPlot
# using SparseArrays
# using Random
# using LinearAlgebra
# using Statistics
# using Arpack

m=2
n=4
gr = grid2(m,n)
(x,y) = grid2coords(4,3)
pyplot(); plot_graph(gr,x,y;dots=false)

Matrix(lap(gr))

Matrix(lap(grid2(m,n)))


cgr = complete_graph(m*n)
Matrix(lap(cgr))

Matrix(lap(grid2(m,n;isotropy=.2)))

Matrix(lap(path_graph(m*n)))
