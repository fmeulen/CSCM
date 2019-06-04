"""
Conversion vector and matrix format
Example usage
    a =  reshape(1:6,2,3)
    vec2mat(mat2vec(a),2,3)-a
"""
mat2vec(x) = vec(x) # colunmwise filling

vec2mat(x,m,n) = reshape(x,m,n)

sigmoid(x) = exp(x)/(1+exp(x))

invlogit(x) = sigmoid.(x)/sum(sigmoid.(x))
#invlogit(x) = exp.(x)/sum(exp.(x))

extractind(w,i) = map(x->x[i],w)

indbin(x_,binx) = (x_>=binx[end]) ? error("bin does not exist") : findfirst(x -> x > x_, binx) -1  # look up in which bin x_ is located

graphlaplacian(m,n) = Matrix(lap(grid2(m,n))) + I/(m*n)^2


mutable struct censoringinfo
    fracarea  # keep track of fraction of bin areas
    ind   # corresponding indices
end


function construct_censoringinfo(t,y,binx,biny,z,ind_yknown)
    m = length(binx) - 1
    n = length(biny) - 1
    N = length(t)
    z = zeros(Int64,N)
    z[ind_yknown] .= 1 # so if y is known z=1
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
    ci
end

bernpar(θ,cik) = dot(θ[cik.ind], cik.fracarea)

@model GraphLaplacianModel(z,ci,L) = begin
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormalCanon(L*τ)
    θ = invlogit(H)
    for k in eachindex(z)
        z[k] ~ Bernoulli(bernpar(θ,ci[k]))
    end
end

@model DirichletFullModel(z,L) = begin
    τ ~ InverseGamma(.1,.1)
    H ~ MvNormalCanon(L*τ)
    θ = invlogit(H)
    for k in ind
        z[k] ~ Bernoulli(θ[k])
    end
end

function graphlaplacian_inference(ci, L, gg)
    model_gl = GraphLaplacianModel(ones(Int8,length(ci)),ci,L)
    Turing.sample(model_gl, gg)  # here also possible to use Gibbs sampling where tau and H are iteratively updated.
end

function dirichletfull_inference(counts,θ,L,gg)
    model_gl = DirichletFullModel(counts,θ,L)
    Turing.sample(model_gl, gg)  # here also possible to use Gibbs sampling where tau and H are iteratively updated.
end


"""
Compute inverse graph Laplacian, with power parameter (ρ) fixed to one.
"""
function invgraphlaplacian(m,n,τ) # order columnwise
    kol1 = [2; fill(3,n-2); 2]
    diagD = vcat(kol1, repeat(kol1 .+ 1,m-2),kol1)
    k = m * n
    L = diagm(0=>diagD .+ 1/k^2, 1=> fill(-1,k-1), -1 => fill(-1,k-1))   # graph Laplacian matrix
    Σ = τ^(-1) * inv(L)
    (Σ+Σ')/2
end

"""
Generate data for current status continuous mark model.
"""
function gendata(truedatagen, N, θ)  # θ is par for GaussianCopula copula
    t = sqrt.(rand(N))
    if truedatagen=="x+y"
        x = -0.5 .+ 0.5 *sqrt.(1 .+ 8*rand(N))
        y = -x .+ sqrt.(x.^2 .+ (2x .+ 1) .* rand(N))
    end
    if truedatagen=="uniform"
        x = rand(N)
        y = rand(N)
    end
    if truedatagen=="(3/8)(x2+y)"
        u = rand(N)
        x = (2u .+ sqrt.(4*u.^2 .+ 1)).^(1/3) .- abs.(2u .- sqrt.(4*u.^2 .+ 1)).^(1/3)
        y = -x.^2 + sqrt.(x.^4 .+ 4*(x.^2 .+ 1) .* rand(N))
    end
    if truedatagen=="GaussianCopula"
        out = rand(GaussianCopula(θ),N)
        x = out[:,1]
        y = out[:,2]
    end
    # find indices of case where x<t (y observed) and x>= t (y not observed)
    ind_yknown = findall(x.<t)
    ind_yunknown = findall(x.>=t)
    # So as available data we are given ind_yknown, ind_yunknown, t and y[ind_yknown]

    x, y, t, ind_yknown, ind_yunknown
end

function Ftrue(x,y,truedatagen, θ)  # θ is par for GaussianCopula copula
    if truedatagen=="uniform"
        out = x+y
    end
    if truedatagen=="x+y"
        out = 0.5*x^2*y + 0.5*x*y^2
    end
    if truedatagen=="(3/8)(x2+y)"
        out = (3/8)*(y*(x^3)/3 + 0.5*x*y^2)
    end
    # if truedatagen=="GaussianCopula"
    #     if abs(x*y) <10^(-2)
    #         out = 0.0
    #     else
    #         out = max((x^(-θ)+y^(-θ)-1)^(-1/θ),0)
    #     end
    # end
    if truedatagen=="GaussianCopula"
        th = θ
        @rput x
        @rput y
        @rput th
        R"""
            library(copula)
            cop = normalCopula(th,2)
            out = pCopula(c(x, y),cop)
        """
        @rget out
    end
    out
end

function ftrue(x,y,truedatagen, θ)  # θ is par for GaussianCopula copula
    if truedatagen=="uniform"
        out = 1
    end
    if truedatagen=="x+y"
        out = x+y
    end
    if truedatagen=="(3/8)(x2+y)"
        out = (3/8)*(x^2 + y)
    end
    # if truedatagen=="GaussianCopula"
    #     if abs(x*y) <10^(-2)
    #         out = 0.0
    #     else
    #         #out = (θ+1) * (x^(-θ) + y^(-θ) -1)^(-1/θ-2)  * (x*y)^(-θ-1)
    #         logout = log(θ+1) - (θ^(-1) + 2)* log(x^(-θ) + y^(-θ) -1) - (θ+1)* (log(x) + log(y))
    #         out = exp(logout)
    #     end
    if truedatagen=="GaussianCopula"
        th = θ
        @rput x
        @rput y
        @rput th
        R"""
            library(copula)
            cop = normalCopula(th,2)
            out = dCopula(c(x, y),cop)
        """
        @rget out
    end
    out
end

"""
Compute the true bin probabilities for data generated by 'truedatagen'.
Probabilities are stacked into a vector where the inner loop goes over the y-direction.
"""
function binprobtrue(binx,biny,truedatagen,θ)
    m,n = length(binx)-1, length(biny)-1
    out = Float64[]
    for i in 1:m
        for j in 1:n
            val = Ftrue(binx[i+1],biny[j+1], truedatagen,θ) -
                        Ftrue(binx[i+1],biny[j], truedatagen,θ) -
                        Ftrue(binx[i],biny[j+1], truedatagen,θ) +
                        Ftrue(binx[i],biny[j], truedatagen,θ)
            push!(out, val)
        end
    end
    out
end

function denstrue(x,y,truedatagen,θ)
    out = Float64[]
    for i in eachindex(x)
        for j in eachindex(y)
            val = ftrue(x[i],y[j], truedatagen,θ)
            push!(out, val)
        end
    end
    out
end

function f_piecewise_constant(x,y,binx,biny,dweights)
    n = length(binx)-1
    ix = indbin(x,binx)
    iy = indbin(y,biny)
    dweights[(ix-1)*n + iy]
end

function dens_piecewise_constant(gridx,gridy,binx,biny,dweights)
    out = Float64[]
    for i in eachindex(gridx)
        for j in eachindex(gridy)
            val = f_piecewise_constant(gridx[i],gridy[j],binx,biny,dweights)
            push!(out, val)
        end
    end
    out
end

function plot_truedensity(d)
    @rput d
    R"""
        library(ggplot2)
        library(tidyverse)
        ggplot(data=d,aes(x, y,z=dtrue)) +
            geom_tile(aes(fill=dtrue))+
            #stat_contour(bins=6,aes(x,y,z=dtrue), color="black", size=0.6)+
            scale_fill_gradient2(low="white", high="black")#+geom_contour(binwidth = 0.5)
            ggtitle("Data generating density") +
             theme_light()

        ggsave("./out/truedensitydensity.pdf",width=4, height=4)

        # geom_contour(binwidth = 0.1) +
        #     #scale_fill_gradient2(limits=c(mincol_lim, maxcol_lim))+
        #     #scale_fill_gradient2()+
    """
end


# Example with data from GaussianCopula copula

abstract type Copula
end

struct GaussianCopula <: Copula
    θ
end


"""
 rand(fam::GaussianCopula,n::Integer)

 simulate bivariate data from GaussianCopula family (returns 2 x n matrix so that each column contains an independent realisation)

 x=rand(GaussianCopula(2.0),5000)
"""
function Base.rand(fam::GaussianCopula,n::Integer)
    th = fam.θ
    @rput th
    @rput n
    R"""
        library(copula)
        cop=normalCopula(th,2)
        out = rCopula(n,copula = cop)
    """
    @rget out
    out
    # θ = fam.θ
    # res = zeros(2,n)
    # for i in 1:n
    #     U = rand()
    #     S = U^(-θ) - 1
    #     res[:,i] = [U, ((1+S) * rand()^(-θ/(1+θ)) - S)^(-1/θ)]
    # end
    # res
end




function prep_plotting_p(pest_dir,pest_gl,truedatagen,binx,biny,titel,θ)
    m = length(binx)-1 ; n = length(biny)-1
    xx = repeat(binx[2:end],inner=n)
    yy = repeat(biny[2:end],outer=m)
    ptrue = binprobtrue(binx,biny,truedatagen,θ)  # true bin probabilities
    d = DataFrame(pest=vcat(pest_dir,pest_gl), ptrue=repeat(ptrue,2),
        x=repeat(xx,2), y=repeat(yy,2),
        type = repeat(["Dirichlet", "Graph Laplacian"],inner=m*n))
    CSV.write("./out/"*titel*"binprob.csv",d)
    d
end

function plotting_p(d)
    @rput d
    R"""
        library(ggplot2)
        library(tidyverse)
        d <- d %>% mutate(error=pest-ptrue)
        ggplot(data=d,aes(x, y, fill=error)) + geom_tile() +
        facet_wrap(~type) +
        scale_fill_gradient2()+xlab("")+ylab("")+
                 theme_light()+ coord_fixed()
        ggsave("./out/binprobs.pdf",width=8.2, height=4)
        errs = d %>% group_by(type) %>%   summarise(error=sqrt(sum(error^2)))
    """
    @rget errs
    errs
end

"""
Compare estimated piecewise constant probability density function, specified via
dweights, to true density computed using treudatagen
"""
function prep_plotting_d(dweights_dir,dweights_gl,truedatagen,binx,biny,θ ;gridN=200)
    # Asses error by binning for true pdf
    minx, maxx = extrema(binx)
    miny, maxy = extrema(biny)
    gridx = range(minx,stop=maxx-0.001,length=gridN)
    gridy = range(miny,stop=maxy-0.001,length=gridN)

    d_dir = dens_piecewise_constant(gridx,gridy,binx,biny,dweights_dir)
    d_gl = dens_piecewise_constant(gridx,gridy,binx,biny,dweights_gl)
    dtrue = denstrue(gridx,gridy,truedatagen,θ)
    x =repeat(gridx,inner=gridN)
    y=repeat(gridy,outer=gridN)
    d = DataFrame(dest=vcat(d_dir,d_gl),dtrue=repeat(dtrue,2),
            x=repeat(x,2), y=repeat(y,2),
            type = repeat(["Dirichlet", "Graph Laplacian"],inner=length(d_dir)))

    CSV.write("./out/density.csv",d)
    d
end

function plotting_d(d)
    @rput d
    R"""
        library(ggplot2)
        library(tidyverse)
        d <- d %>% mutate(error=dest-dtrue)  #error is estimted - true
        ggplot(data=d,aes(x, y, fill=error)) + geom_tile() +
        facet_wrap(~type) +
        scale_fill_gradient2()+xlab("")+ylab("")+
                 theme_light()+ coord_fixed()
        ggsave("./out/dens.pdf",width=8.2, height=4)
        errs = d %>% group_by(type) %>%   summarise(error=sqrt(sum(error^2)))
    """
    @rget errs
    errs
end
