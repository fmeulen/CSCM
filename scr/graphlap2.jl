struct censoringinfo{S<:Number, T<:Number}
    fracarea::Vector{S}         # keep track of fraction of bin areas
    ind::Vector{T}              # corresponding indices
end

graphlaplacian(m,n) = Matrix(lap(grid2(m,n))) + I/(m*n)^2


"""
    mapsimplex(x)

apply logistic funtion to each element of vector `x` and scale such that the
elements sum to 1
"""
function mapsimplex(x)
    y = logistic.(x)
    y/sum(y)
end

@model GraphLaplacianMod(ci,L) = begin
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormalCanon(L*τ)
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
    t: observed times
    ind_yknown: vector of indices that correspond to those times in t where y is observed
    ind_yunknown: vector of indices that correspond to those times in t where y is unobserved
    y: vector of observed marks (elements corresponding to times where the mark y is not observed can be specified arbitrarily,
    for example z zero)
    binx: bin grid in x direction
    biny: bin grid in y direction
    alg: algorithm for sampling

    ci::censoringinfo, chn::Chains
"""
function sample_graphlap(t,ind_yknown, ind_yunknown, y, (binx, biny), ITER; alg=HMC(0.1, 5))
    nsample = length(t)
    m = length(binx) - 1
    n = length(biny) - 1
    # construct censoringinfo
    ci = Vector{censoringinfo}(undef,nsample)
    for k ∈ ind_yknown
        it = indbin(t[k],binx)
        iy = indbin(y[k],biny)
        fa =  [ (min(binx[i+1],t[k])-binx[i])/(binx[i+1]-binx[i]) for i ∈ 1:it]
        ind = [iy + ℓ*n for ℓ ∈ 0:(it-1)]
        ci[k] = censoringinfo(fa, ind)
    end
    for k ∈ ind_yunknown
        it = indbin(t[k],binx)
        fa = [(binx[i+1]-max(t[k],binx[i]))/(binx[i+1] - binx[i])  for i ∈ it:m for j ∈ 1:n]
        ind = collect(((it-1)*n+1):(m*n))
        ci[k] = censoringinfo(fa,ind)
    end
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

    p1 = plot(out1, label = "H[1]")
    p2 = plot(out2, label = "H[2]")
    p10 = plot(out10, label="H[10]")
    pτ = plot(outτ, label="τ")
    l = @layout [a b; c d]
    Plots.plot(p1, p2, p10, pτ, layout=l)
end
