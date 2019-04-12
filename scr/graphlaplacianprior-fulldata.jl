countsfull = zeros(Int64,m,n)
for i in 1:N
    ix = indbin(x[i], binx)
    iy = indbin(y[i], biny)
    countsfull[ix,iy] += 1
end
fdata = vec(countsfull) # full data

@model GraphLaplacionModel(z,N,Linv1) = begin
    k = size(Linv1)[1]
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormal(zeros(k),Linv1/τ)
    θ = invlogit(H)
    z ~ Multinomial(N,θ)
end

Linv1 = invgraphlaplacian(m,n,1)
model_glc = GraphLaplacionModel(fdata,N,Linv1)
chn_glc = sample(model_glc, gg)

iterates_glc = chn_glc.value[BI:end,:,1]
# Hpostmean_glc = mean(iterates_glc,dims=1)[1:m*n]  # rowmeans
# θpostmean_glc = invlogit(Hpostmean_glc)
Hpostmean_glc = chn_glc.value[BI:end,1:m*n,1]
θiterates_glc = mapslices(invlogit, Hpostmean_glc,dims=2)
θpostmean_glc = vec(mean(θiterates_glc,dims=1))


plotly(); plot(chn_glc)
describe(chn_glc)
plotting(θpostmean,truedatagen,binarea,binx,biny,"fulldata graph Laplacian","contour fulldata graph Laplacian")
