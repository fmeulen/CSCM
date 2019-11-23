#cd("/Users/Frank/Sync/DOCUMENTS/onderzoek/code/currentstatuscontinuousmarks")
using Distributions, Test, Statistics, Random, LinearAlgebra
using DelimitedFiles,  DataFrames, RCall
using Turing
using StatsPlots
using DynamicHMC
using CSV, Statistics
using Laplacians
#@rlibrary copula

#----------------------------------------------------------------------------------------------
workdir = @__DIR__
println(workdir)
cd(workdir)
include("funcdefs.jl")
Random.seed!(1234)

#----------------------------------------------------------------------------------------------
# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)","GaussianCopula"][4]
NSAMPLE = 500 # sample size
θcopula = -0.65 # par for GaussianCopula copula
x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen, NSAMPLE, θcopula)

#----------------------------------------------------------------------------------------------
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

#----------------------------------------------------------------------------------------------
ITER_DIR = 1000 #nr of iterations
BI_DIR = 100 #nr of burnin iterations

ITER_GL = 10
BI_GL = div(ITER_GL,3)

# Dirichlet prior
ps = 0.1
@time iterates_dir, θpostmean_dir = sample_dir(t,ind_yknown, y,binx,biny, BI_DIR, ITER_DIR; priorscale = ps)

# Graph Laplacian prior
samplers=[HMC(0.1, 5) HMC(0.1, 2)  DynamicNUTS()]
sp = HMC(0.1, 5)#samplers[1]

@time iterates_gl, Hiterates_gl, θiterates_gl, θpostmean_gl = sample_graphlap(t,ind_yknown, y,binx,biny, ITER_GL, BI_GL; sp=sp)

# ensure that there is a directory called "out" in the working directory
make_traceplots = false
include("write_info.jl")
