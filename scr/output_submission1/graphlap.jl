mutable struct censoringinfo
    fracarea  # keep track of fraction of bin areas
    ind   # corresponding indices
end

# Turing hierarchical model
@model GraphLaplacianModel(z,ci,L) = begin
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormalCanon(L*τ)
    θ = invlogit(H)
    for k in eachindex(z)
        z[k] ~ Bernoulli(sum(θ[ci[k].ind].* ci[k].fracarea))
    end
end

"""
    t: observed times
    ind_yknown: vector of indices that correspond to those times in t where y is observed
    y: vector of observed marks (elements corresponding to times where the mark y is not observed can be specified arbitrarily,
    for example z zero)
    binx: bin grid in x direction
    biny: bin grid in y direction
    BI: number of samples considered burnIn
    sampler: sampler

    Returns:     iterates, Hiterates, θiterates, θpostmean
    Here 'iterates' contains all iterates, whereas for
    Hiterates, θiterates, θpostmean burnin samples have been removed
"""
function sample_graphlap(t,ind_yknown, y,binx,biny, ITER_GL, BI_GL; sampler=HMC(0.1, 5))
    binarea = (binx[2]-binx[1]) * (biny[2]-biny[1]) # the same for all bins
    NSAMPLE = length(t)
    zz = zeros(Int64,NSAMPLE)
    zz[ind_yknown] .= 1 # so if y is known z=1
    m = length(binx) - 1
    n = length(biny) - 1
    ci = Vector{censoringinfo}(undef,NSAMPLE)

    for k in 1:NSAMPLE
        it = indbin(t[k],binx)
        if zz[k]==1
            iy = indbin(y[k],biny)
            area = [ min(binx[i+1],t[k])-binx[i] for i in 1:it] * (biny[iy+1] - biny[iy])
            ind = [iy + ℓ*n for ℓ in 0:(it-1)]
        else
            area = [(binx[i+1]-max(t[k],binx[i])) * (biny[j+1] - biny[j])  for i in it:m for j in 1:n]
            ind = ((it-1)*n+1):(m*n)
        end
        ci[k] = censoringinfo(area/binarea,ind)
    end

    L = graphlaplacian(m,n) # graph Laplacian with τ=1
    model = GraphLaplacianModel(ones(Int8,NSAMPLE),ci,L)
    chn = Turing.sample(model, sampler,ITER_GL)
    # here also possible to use Gibbs sampling where tau and H are iteratively updated.

    iterates = chn.value[1:end,:,1]
    Hiterates = chn.value[BI_GL:end,1:m*n,1]
    θiterates = mapslices(invlogit, Hiterates,dims=2) # transform iterates back to θ
    θpostmean = vec(mean(θiterates,dims=1))
    iterates, Hiterates, θiterates, θpostmean
end
