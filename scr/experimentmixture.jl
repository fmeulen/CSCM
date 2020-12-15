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

include("funcdefs.jl")
include("mcmc.jl")

Random.seed!(1234)

dat = load("data_Mixture.jld2")
ind_yunknown, t, nsample, x, dist, ind_yknown, y =
                dat["ind_yunknown"], dat["t"], dat["nsample"], dat["x"], dat["dist"], dat["ind_yknown"], dat["y"]

mkpath("./out/mixture/mnsmall")
mkpath("./out/mixture/mnmedium")
mkpath("./out/mixture/mnlarge")



IT = 20_000
BI = div(IT,3)
bi = BI:IT

# Prior on τ
Πdirτ = Exponential(1.0)
Πτ = Exponential(1.0)
outdir_options = ["./out/mixture/mnsmall", "./out/mixture/mnmedium", "./out/mixture/mnlarge"]
bins_options = [Bins(dist, 5, 10), Bins(dist, 25, 50), Bins(dist, 100, 200)]

for i ∈ 1:3
    bins = bins_options[i]
    ci = construct_censoringinfo(t, y, ind_yknown, ind_yunknown, bins)
    # Dirichlet prior
    @time θdir, τdir, accdir = dirichlet(ci, bins, IT, Πdirτ)
    @show accdir/IT

    # pCN for Laplacian prior
    @time θgl, τgl, accgl, ρ = pcn(ci, bins, IT, Πτ; ρ=.96, δ=0.6)
    @show accgl/IT

    outdir = outdir_options[i]
    processoutput(θdir, τdir, accdir, θgl, τgl, accgl, ρ,
            dist, nsample, bins, BI, IT, outdir)
end
