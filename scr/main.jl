using Distributions, Random, LinearAlgebra
using DataFrames, RCall
using StatsFuns
using CSV
using Laplacians
using Plots
using PDMats
using FileIO


wd = @__DIR__
cd(wd)
include("funcdefs.jl")
include("mcmc.jl")

mkpath("./out")
Random.seed!(1234)

# Read data
dat = load("data1.jld2")
ind_yunknown, t, nsample, x, dist, ind_yknown, y = dat["ind_yunknown"], dat["t"], dat["nsample"], dat["x"], dat["dist"], dat["ind_yknown"], dat["y"]

# Compute bins
bins = Bins(dist, 10, 20)

# Combine observations to type CensoringInfo (note that y[ind_yunknown] can be anything)
ci = construct_censoringinfo(t, y, ind_yknown, ind_yunknown, bins)

IT = 10_000 # nr of iterations for Dirichlet prior
BI = div(IT,3) # nr of burnin iters
bi = BI:IT

# Prior on τ
Πτ = Exponential(1.0)

# Dirichlet prior
@time θdir, τdir, accdir = dirichlet(ci, bins, IT, Πτ)
@show accdir/IT

# Pcn for Laplacian prior
@time θgl, τgl, accgl, ρ = pcn(ci, bins, IT, Πτ; ρ=.96, δ=0.6)
@show accgl/IT

include("write_info.jl")
