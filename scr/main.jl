cd("/Users/Frank/Sync/DOCUMENTS/onderzoek/code/currentstatuscontinuousmarks")
using Distributions, Test, Statistics, Random, LinearAlgebra
using DelimitedFiles,  DataFrames, RCall
using Turing, StatsPlots, DynamicHMC, CSV, Statistics
using Laplacians

#@rlibrary copula

include("funcdefs.jl")

Random.seed!(1234)

# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)","GaussianCopula"][4]
N = 500 # sample size
θcopula = -0.65 # par for GaussianCopula copula
x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen, N, θcopula)

# for Dirichlet prior specify hyperparameter
priorscale = 0.1

# Set bins (m and n are nr of bins in horizontal and vertical directions respectively)
minx = 0.0; miny = 0.0; maxx = 1.0
m = 10

if truedatagen=="(3/8)(x2+y)"
    n = m
    maxy = 2.0
else
    n = m
    maxy = 1.0
end
binx = range(minx, stop=maxx, length=m+1)
biny = range(miny, stop=maxy, length=n+1)
binarea = (binx[2]-binx[1]) * (biny[2]-biny[1]) # the same for all bins


ITER = 1000 #nr of iterations
BI = 100 #nr of burnin iterations

include("dirichletprior.jl")

samplers=[HMC(100, 0.1, 5) HMC(ITER, 0.1, 2)  DynamicNUTS(100)  PG(10,1000) SGLD(1000, 0.5) SGHMC(100,0.01,0.01)]
gg = samplers[1]
#include("graphlaplacianprior.jl")
include("graphlaplacianprior.jl")

# ensure that there is a directory called "out" in the working directory
include("write_info.jl")
