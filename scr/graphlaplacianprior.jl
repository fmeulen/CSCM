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

@model GraphLaplacianModel(z,ci,L) = begin
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormalCanon(L*τ)
    θ = invlogit(H)
    for k in eachindex(z)
        z[k] ~ Bernoulli(sum(θ[ci[k].ind].* ci[k].fracarea))
    end
end

L = graphlaplacian(m,n) # graph Laplacian with τ=1
model_gl = GraphLaplacianModel(ones(Int8,N),ci,L)
chn_gl = Turing.sample(model_gl, gg)
# here also possible to use Gibbs sampling where tau and H are iteratively updated.

iterates_gl = chn_gl.value[1:end,:,1]
BI = 1
Hpostmean_gl = chn_gl.value[BI:end,1:m*n,1]
θiterates_gl = mapslices(invlogit, Hpostmean_gl,dims=2) # transform iterates back to θ
θpostmean_gl = vec(mean(θiterates_gl,dims=1))

plotly(); plot(chn_gl)
describe(chn_gl)

trace_df = DataFrame(th1=θiterates_gl[:,1],th2=θiterates_gl[:,2],th10=θiterates_gl[:,10],
th11=θiterates_gl[:,11])
nr, nc = size(trace_df)
@rput trace_df
@rput nr; @rput nc
R"""
library(tidyverse)
trace_df %>% gather(value=value, key=coefficient) %>% mutate(iterate=rep(1:nr,nc)) %>%
ggplot(aes(x=iterate,y=value)) + geom_path() + facet_wrap(~coefficient,scales='free')
"""
