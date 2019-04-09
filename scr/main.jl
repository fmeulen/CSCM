cd("/Users/Frank/Sync/DOCUMENTS/onderzoek/code/currentstatuscontinuousmarks")
using Distributions, Test, Statistics, Random, LinearAlgebra
using DelimitedFiles,  DataFrames, RCall
using Turing,  StatsPlots, DynamicHMC, CSV

include("funcdefs.jl")

Random.seed!(1234)

# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)"][3]
N = 1000 # sample size
x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen, N)

# for Dirichlet prior specify hyperparameter
priorscale = 0.1

# Set bins (m and n are nr of bins in horizontal and vertical directions respectively)
minx = 0.0; miny = 0.0; maxx = 1.0
m = 5

if truedatagen=="(3/8)(x2+y)"
    n = 2m
    maxy = 2.0
else
    n = m
    maxy = 1.0
end
binx = range(minx, stop=maxx, length=m+1)
biny = range(miny, stop=maxy, length=n+1)
binarea = (binx[2]-binx[1]) * (biny[2]-biny[1]) # the same for all bins


ITER = 2000 #nr of iterations
BI = 200 #nr of burnin iterations

include("dirichletprior.jl")

samplers=[HMC(ITER, 0.1, 10)  DynamicNUTS(100)  PG(10,1000)]
gg = samplers[2]
#include("graphlaplacianprior.jl")
include("graphlaplacianprior_po.jl")

# ensure that there is a directory called "out" in the working directory
include("write_info.jl")
