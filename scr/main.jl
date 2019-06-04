#cd("/Users/Frank/Sync/DOCUMENTS/onderzoek/code/currentstatuscontinuousmarks")
cd("/Users/Frank/github/CSCM/scr")
using Distributions, Test, Statistics, Random, LinearAlgebra
using DelimitedFiles,  DataFrames, RCall
using Turing, StatsPlots, DynamicHMC, CSV, Statistics
using Laplacians

#@rlibrary copula

include("funcdefs.jl")

Random.seed!(1234)

# Sample data (available data consist of (ind_yknown, ind_yunknown, t, y[ind_yknown])
truedatagen =["uniform","x+y","(3/8)(x2+y)","GaussianCopula"][2]
N = 10 # sample size
θcopula = -0.7 # par for GaussianCopula copula
x, y, t, ind_yknown, ind_yunknown = gendata(truedatagen, N, θcopula)

# for Dirichlet prior specify hyperparameter
priorscale = 0.1

# Set bins (m and n are nr of bins in horizontal and vertical directions respectively)
minx = 0.0; miny = 0.0; maxx = 1.0
m = 5
if truedatagen=="(3/8)(x2+y)"
    n = m
    maxy = 2.0
else
    n = m
    maxy = 1.0
end
binx = range(minx, stop=maxx, length=m+1)
biny = range(miny, stop=maxy, length=n+1)
binarea = (binx[2]-binx[1]) * (biny[2]-biny[1])  # the same for all bins


ITER = 5000  #nr of iterations for Dirichlet prior
BI = div(ITER,2)#nr of burnin iterations

############### Dirichlet prior ###########################
include("dirichletprior.jl")

############### graph Laplacian prior ###########################

#samplers=[HMC(100, 0.05, 5) HMC(ITER, 0.1, 2)  DynamicNUTS(100)  PG(10,1000) SGLD(1000, 0.5) SGHMC(100,0.01,0.01)]
#gg = samplers[5]

ci = construct_censoringinfo(t,y,binx,biny,z, ind_yknown)
L = graphlaplacian(m,n) # graph Laplacian with τ=1

Turing.setadbackend(:forward_diff) #Turing.setadbackend(:reverse_diff)

ITER_NUTS = 2000
BI_NUTS = 100
gg = DynamicNUTS(ITER_NUTS)

offs = 1/(m*n)^2
@time chn = graphlaplacian_inference(ci, L + offs*I, gg)

## summarise output
ind = BI_NUTS:ITER_NUTS
itgl = Float64.(chn.value[:,:,1])
τgl = Float64.(itgl[ind,end])
Hgl = itgl[ind,1:m*n] # posterior mean
θitgl = mapslices(invlogit, Hgl,dims=2) # iterates (transform iterates back to θ)
θpm_gl = vec(mean(θitgl,dims=1))
τpm_gl = mean(τgl)

plotly(); plot(chn)
#describe(chn_gl)

# ensure that there is a directory called "out" in the working directory
include("write_info.jl")
