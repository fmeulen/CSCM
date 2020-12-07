using Distributions, Random, LinearAlgebra
#using Test, Statistics
using DelimitedFiles,  DataFrames, RCall
using Turing
using StatsFuns
using DynamicHMC
using CSV, Statistics
using Laplacians
using Setfield
using MCMCChains
#using ReverseDiff
using Zygote
using Plots
using Optim

@rlibrary(transport)

workdir = @__DIR__
cd(workdir)
include("funcdefs.jl")
include("dirichlet.jl")
include("graphlap2.jl")
include("pcn.jl")
Random.seed!(1234)

#Turing.setadbackend(:reversediff)
Turing.setadbackend(:zygote)

# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)","multimodal"][4]
nsample = 1000 # sample size
x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen,nsample)

# Set bins (m and n are nr of bins in horizontal and vertical directions respectively)
minx = 0.0; miny = 0.0; maxx = 1.0
m = 35
if truedatagen=="(3/8)(x2+y)"
    n = m
    maxy = 2.0
else
    n = m
    maxy = 1.0
end
binx = range(minx, stop=maxx, length=m+1)
biny = range(miny, stop=maxy, length=n+1)


IT = 10_0 # nr of iterations for Dirichlet prior
BI = div(IT,3) # nr of burnin iters
bi = BI:IT

# Dirichlet prior
ps = 0.1
@time θdir = sample_dir(t,ind_yknown, y, (binx, biny), IT; priorscale = ps)

# Graph Laplacian prior
if false
    ITERgl = 5_000 # nr of iterations for GraphLaplacian prior
    BIgl = div(ITERgl,3) # nr of burnin iters
    bi_gl = BIgl:ITERgl
    samplers=[HMC(0.1, 5), HMC(0.2, 20), HMC(0.05,20),  DynamicNUTS()]
    sp = samplers[2]
    @time  ci, chn, τ, H, θgl=  sample_graphlap(t,ind_yknown, ind_yunknown, y, (binx, biny), ITERgl; alg=sp)
end

@time θsave, τsave, acc = pcn(t,ind_yknown, y, (binx, biny), IT)


# ensure that there is a directory called "out" in the working directory
include("write_info.jl")
