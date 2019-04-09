z = zeros(Int64,N)
z[ind_yknown] .= 1 # so if y is known z=1

mutable struct censoringinfo
    fracarea  # keep track of fraction of bin areas
    ind   # corresponding indices
end

ci = Vector{censoringinfo}(undef,N)

for k in 1:N
    it = indbin(t[k],binx)
    if z[k]==1
        iy = indbin(y[k],biny)
        area = [ min(binx[i+1],t[k])-binx[i] for i in 1:it] * (biny[iy+1] - biny[iy])
        ind = [iy + ℓ*n for ℓ in 0:(it-1)]
    else
        area = [(binx[i+1]-max(t[k],binx[i])) * (biny[j+1] - biny[j])  for i in it:m for j in 1:n]
        ind = ((it-1)*n+1):(m*n)
    end
    ci[k] = censoringinfo(area/binarea,ind)
end

@model GraphLaplacionModel_po(z,N,ci,Linv1) = begin
    k = size(Linv1)[1]
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormal(zeros(k),Linv1/τ)
    θ = invlogit(H)
    for k in 1:N
            z[k] ~ Bernoulli(sum(θ[ci[k].ind].* ci[k].fracarea))
    end
end

Linv1 = invgraphlaplacian(m,n,1) # graph Laplacian with τ=1
model_gl = GraphLaplacionModel_po(ones(Int8,N),N,ci,Linv1)
chn_gl = sample(model_gl, gg)


# here also possible to use Gibbs sampling where tau and H are iteratively updated.


iterates_gl = chn_gl.value[1:end,:,1]
# Hpostmean_gl = mean(iterates_gl,dims=1)[1:m*n]  # rowmeans
# θpostmean_gl = invlogit(Hpostmean_gl)
BI = 1
Hpostmean_gl = chn_gl.value[BI:end,1:m*n,1]
θiterates_gl = mapslices(invlogit, Hpostmean_gl,dims=2) # transform iterates back to θ
θpostmean_gl = vec(mean(θiterates_gl,dims=1))

plotly(); plot(chn_gl)
describe(chn_gl)
plotting(θpostmean_gl,truedatagen,binarea,binx,biny,"graphLaplacian","contourGraphLaplacian")
