# Turing implementation

using Turing
Turing.setadbackend(:zygote)
using DynamicHMC
using MCMCChains
using Zygote


@model GraphLaplacianMod(ci,L) = begin
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormalCanon(L/τ)
    Turing.@addlogprob! loglik(H, ci)
end

function loglik(H, ci)
    θ = softmax(H)
    ll = 0.0
    @inbounds for i ∈ eachindex(ci)
        γ = ci[i]
        ll += log(dot(θ[γ.ind], γ.fracarea))
    end
    ll
end

"""
    sample_graphlap(ci, (m, n), IT; alg=HMC(0.1, 5))

ci:: CensoringInfo
(m, n): number of horizontal and vertical bins
IT: number of iterations

Returns
chn::Chains
τ: iterates
H: iterates
θ: iterates
"""
function sample_graphlap(ci, (m, n), IT; alg=HMC(0.1, 5))
    ci = construct_censoringinfo(t, (binx,biny), ind_yknown, ind_yunknown)
    # define model
    L = graphlaplacian(m,n) # graph Laplacian with τ=1
    model = GraphLaplacianMod(ci,L)
    # run sampler
    chn = Turing.sample(model, alg, IT)

    # extracting H, θ, τ
    τ = chn[:,:τ,1].data
    H = chn.value[:,1:m*n,1]
    θ = mapslices(softmax, H,dims=2)

    chn, τ, H, θ
end

function traceplots(chn)
    outτ = chn[:,:τ,:]
    out1 =  chn[:,Symbol("H[1]"),:]
    out2 =  chn[:,Symbol("H[2]"),:]
    out10 = chn[:,Symbol("H[10]"),:]

    p1 = Plots.plot(out1, label = "H[1]")
    p2 = Plots.plot(out2, label = "H[2]")
    p10 = Plots.plot(out10, label="H[10]")
    pτ = Plots.plot(outτ, label="τ")
    l = @layout [a b; c d]
    Plots.plot(p1, p2, p10, pτ, layout=l)
end


## functions calls, actual programme

ITERgl = 5_000 # nr of iterations for GraphLaplacian prior
BIgl = div(ITERgl,3) # nr of burnin iters
bi_gl = BIgl:ITERgl
samplers=[HMC(0.1, 5), HMC(0.2, 20), HMC(0.05,20),  DynamicNUTS()]
sp = samplers[2]
@time  chn, τ, H, θgl=  sample_graphlap(ci, (m,n), ITERgl; alg=sp)


θ̄gl = vec(mean(θgl[bi_gl,:], dims=1))  # posterior mean graphlap using Turing
p = traceplots(chn) # for Turing output
