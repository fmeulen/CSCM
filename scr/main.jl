using Distributions, Random, LinearAlgebra
#using Test, Statistics
using DelimitedFiles,  DataFrames, RCall
using Turing
using StatsPlots
using StatsFuns
using DynamicHMC
using CSV, Statistics
using Laplacians
using Setfield
using MCMCChains

workdir = @__DIR__
println(workdir)
cd(workdir)
include("funcdefs.jl")
include("dirichlet.jl")
include("graphlap2.jl")
Random.seed!(1234)


# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)","GaussianCopula"][3]
nsample = 500 # sample size
θcopula = -0.65 # par for GaussianCopula copula
x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen,nsample, θcopula)

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
#binarea = (binx[2]-binx[1]) * (biny[2]-biny[1]) # the same for all bins




ITERdir = 5000 # nr of iterations for Dirichlet prior
BIdir = 100 # nr of burnin iters
bi_dir = BIdir:ITERdir

ITERgl = 500 # nr of iterations for GraphLaplacian prior
BIgl = div(ITERgl,3) # nr of burnin iters
bi_gl = BIgl:ITERgl

# Dirichlet prior
ps = 0.1
@time θdir = sample_dir(t,ind_yknown, y, (binx, biny), ITERdir; priorscale = ps)

# Graph Laplacian prior
samplers=[HMC(0.1, 5) HMC(0.1, 2)  DynamicNUTS()]
sp = samplers[1]

@time  ci, chn, τ, H, θgl=  sample_graphlap(t,ind_yknown, ind_yunknown, y, (binx, biny), ITERgl; alg=sp)

# ensure that there is a directory called "out" in the working directory
include("write_info.jl")
