using Distributions, Random, LinearAlgebra
#using DelimitedFiles,
using DataFrames, RCall
using StatsFuns
using CSV
#using Statistics
using Laplacians
#using Setfield
using Plots
using PDMats
@rlibrary(transport)


workdir = @__DIR__
cd(workdir)
include("funcdefs.jl")
include("mcmc.jl")

Random.seed!(1234)


# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
#truedatagen =["uniform","x+y","(3/8)(x2+y)","multimodal"][3]
dist = Xplusy()
nsample = 100 # sample size
x, y, t, ind_yknown, ind_yunknown = censdata(dist, nsample)
#x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen,nsample)

# Set bins (m and n are nr of bins in horizontal and vertical directions respectively)

minx = 0.0; miny = 0.0; maxx = 1.0
m = n = 4
maxy = dist==X2plusy() ? 2.0 : 1.0
n = dist==X2plusy() ? 2m : m
binx = range(minx, stop=maxx, length=m+1)
biny = range(miny, stop=maxy, length=n+1)

# combine observations to type CensoringInfo (note that y[ind_yunknown] can be anything)
ci = construct_censoringinfo(t, y, (binx,biny), ind_yknown, ind_yunknown)

IT = 25_00 # nr of iterations for Dirichlet prior
BI = div(IT,3) # nr of burnin iters
bi = BI:IT

# prior on τ
Πτ = Exponential(1.0)

# Dirichlet prior
@time θdir, τdir, accdir = dirichlet(ci, (m, n), IT, Πτ)
@show accdir/IT

# pcn for Laplacian prior
@time θgl, τgl, accgl, ρ = pcn(ci, (m, n), IT, Πτ; ρ=.96, δ=0.6)
@show accgl/IT

# ensure that there is a directory called "out" in the working directory
include("write_info.jl")
