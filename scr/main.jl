using Distributions, Random, LinearAlgebra
using DataFrames, RCall
using StatsFuns
using CSV
using Laplacians
using Plots
using PDMats
using FileIO
using Cubature
using Printf


wd = @__DIR__
cd(wd)
mkpath("./out")
include("funcdefs.jl")
include("mcmc.jl")
include("write_info.jl")

#Random.seed!(1234)

# Simuldate data
dist = Xplusy()
nsample = 100
x, y, t, ind_yknown, ind_yunknown = gencensdata(dist, nsample)

# Compute bins
bins = Bins(dist, 10, 10)

# Combine observations to type CensoringInfo (note that y[ind_yunknown] can be anything)
ci = construct_censoringinfo(t, y, ind_yknown, ind_yunknown, bins)

# Iterations and burnin
IT = 20_000
BI = div(IT,3)


# Dirichlet prior
Πdirτ = Exponential(1.0)   #Uniform(0.01, 1.0)
@time θdir, τdir, accdir = dirichlet(ci, bins, IT, Πdirτ)
@show accdir/IT

# pCN for LNGL-prior
Πτ = Exponential(1.0)
@time θgl, τgl, accgl, ρ = pcn(ci, bins, IT, Πτ; ρ=.96, δ=0.6)
@show accgl/IT

outdir = "./out"
processoutput(θdir, τdir, accdir, θgl, τgl, accgl, ρ,
        dist, nsample, bins, BI, IT, outdir)
