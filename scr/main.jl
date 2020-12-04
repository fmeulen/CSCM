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
@rlibrary copula

workdir = @__DIR__
println(workdir)
cd(workdir)
include("funcdefs.jl")
include("dirichlet.jl")
include("graphlap2.jl")
Random.seed!(1234)

#Turing.setadbackend(:reversediff)
Turing.setadbackend(:zygote)

# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)","GaussianCopula"][3]
nsample = 100 # sample size
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
samplers=[HMC(0.1, 5) HMC(0.1, 10)  DynamicNUTS()]


sp = samplers[2]

@time  ci, chn, τ, H, θgl=  sample_graphlap(t,ind_yknown, ind_yunknown, y, (binx, biny), ITERgl; alg=sp)

# ensure that there is a directory called "out" in the working directory
include("write_info.jl")


 ##FIXME, also compute MLE
optimize(model, MLE(), NelderMead())

# Generate an MAP estimate.
map_estimate = optimize(model, MAP())

# Generate an MAP estimate.
map_estimate = optimize(model, MAP())
# Sample with the MAP estimate as the starting point.
chain = sample(model, NUTS(), 1_000, init_theta = map_estimate.values.array)
