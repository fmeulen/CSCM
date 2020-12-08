mapsimolex(x) = softmax(x)

@model GraphLaplacianMod(ci,L) = begin
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormalCanon(L/τ)
    Turing.@addlogprob! loglik(H, ci)
end

function loglik(H, ci)
    θ = mapsimplex(H)
    ll = 0.0
    @inbounds for i ∈ eachindex(ci)
        γ = ci[i]
        ll += log(dot(θ[γ.ind], γ.fracarea))
    end
    ll
end

"""
    sample_graphlap(t,ind_yknown, ind_yunknown, (binx, biny), ITER; alg=HMC(0.1, 5))

t: observed times
ind_yknown: vector of indices that correspond to those times in t where y is observed
ind_yunknown: vector of indices that correspond to those times in t where y is unobserved
y: vector of observed marks (elements corresponding to times where the mark y is not observed can be specified arbitrarily,
for example z zero)
binx: bin grid in x direction
biny: bin grid in y direction
alg: algorithm for sampling

Returns:
    ci, chn, τ, H, θ

ci::CensoringInfo
chn::Chains
τ: iterates
H: iterates
θ: iterates
"""
function sample_graphlap(t,ind_yknown, ind_yunknown, (binx, biny), ITER; alg=HMC(0.1, 5))
    ci = construct_censoringinfo(t, (binx,biny), ind_yknown, ind_yunknown)
    # define model
    L = graphlaplacian(m,n) # graph Laplacian with τ=1
    model = GraphLaplacianMod(ci,L)
    # run sampler
    chn = Turing.sample(model, alg, ITER)

    # extracting H, θ, τ
    τ = chn[:,:τ,1].data
    H = chn.value[:,1:m*n,1]
    θ = mapslices(mapsimplex, H,dims=2)

    ci, chn, τ, H, θ
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
