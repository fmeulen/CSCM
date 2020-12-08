using Distributions, Random, LinearAlgebra
using DelimitedFiles,  DataFrames, RCall
using StatsFuns
using CSV, Statistics
using Laplacians
using Setfield
using Plots
using PDMats
@rlibrary(transport)

# using Turing
#Turing.setadbackend(:zygote)
# using DynamicHMC
# using MCMCChains
# using Zygote


workdir = @__DIR__
cd(workdir)
include("funcdefs.jl")
include("dirichlet2.jl")
include("pcn.jl")
Random.seed!(1234)



# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)","multimodal"][3]
nsample = 100 # sample size
x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen,nsample)

# Set bins (m and n are nr of bins in horizontal and vertical directions respectively)
minx = 0.0; miny = 0.0; maxx = 1.0
m = n = 40

maxy = truedatagen=="(3/8)(x2+y)" ? 2.0 : 1.0
binx = range(minx, stop=maxx, length=m+1)
biny = range(miny, stop=maxy, length=n+1)


IT = 25_00 # nr of iterations for Dirichlet prior
BI = div(IT,3) # nr of burnin iters
bi = BI:IT

# Dirichlet prior
#@time θdir = sample_dir(t,ind_yknown, y, (binx, biny), IT; priorscale = ps)
@time θdir = dirichlet(t, ind_yknown, y, (binx, biny), IT; priorscale = 0.1)

# Graph Laplacian prior
if false
    include("graphlap2.jl")  # for Turing
    ITERgl = 5_000 # nr of iterations for GraphLaplacian prior
    BIgl = div(ITERgl,3) # nr of burnin iters
    bi_gl = BIgl:ITERgl
    samplers=[HMC(0.1, 5), HMC(0.2, 20), HMC(0.05,20),  DynamicNUTS()]
    sp = samplers[2]
    @time  ci, chn, τ, H, θgl=  sample_graphlap(t,ind_yknown, ind_yunknown, (binx, biny), ITERgl; alg=sp)
end

# pcn for Laplacian prior
@time θsave, τsave, acc, ρ = pcn(t,ind_yknown, y, (binx, biny), IT; ρ=.96, δ=0.4, τinit=100.0)
@show acc/IT

# ensure that there is a directory called "out" in the working directory
include("write_info.jl")
