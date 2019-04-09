mutable struct obs
    t::Float64
    ix::Int64
    iy::Int64
    it::Int64
    obstype::String
    area  # keep track of bin areas or lengths (in case y is known), needed for updating latent data
    indpairs  # only needed for case y is unknown
end

global counts = zeros(Int64,m,n)  # adjust at each iteration
global fulldata = Vector{obs}(undef,N)  # adjust at each iteration
####### Initialise fulldata and counts
for k in 1:N
    it = indbin(t[k],binx)
    if k in ind_yknown
        ix = sample(1:it)
        iy = indbin(y[k],biny)
        obstype="yknown"
        area = [ min(binx[i+1],t[k])-binx[i] for i in 1:it]
        indpairs = 0
    else
        ix = sample(it:m)
        iy = sample(1:n)
        obstype="yunknown"
        area = [(binx[i+1]-max(t[k],binx[i])) * (biny[j+1] - biny[j])  for i in it:m for j in 1:n]
        indpairs =[ [i,j] for i in it:m for j in 1:n]
    end
    counts[ix,iy] +=1
    fulldata[k] = obs(t[k],ix,iy,it,obstype,area,indpairs)
end

priorθ = priorscale *  ones(m,n)

global θ = Vector{Matrix{Float64}}(undef,ITER) # save in each iteration

####### Data augmentation algorithm
for iter in 1:ITER
    # Update weights
    #θ[iter] = reshape( rand(Dirichlet(vec(counts+priorθ))) , m, n)
    θ[iter] = vec2mat( rand(Dirichlet(vec(counts+priorθ))) , m, n)

    # Update latent data
    for k in 1:N#  sample(1:n, div(n,2))#1:n
        ix_ = fulldata[k].ix;  iy_ = fulldata[k].iy; it_ = fulldata[k].it;
        t_ = fulldata[k].t
        obstype_ = fulldata[k].obstype
        area_ = fulldata[k].area

        counts[ix_,iy_] += -1
        if obstype_=="yknown"  # update ix, consider i in 1..it_,  j=iy_
            w = [θ[iter][i,iy_] for i in 1:it_] .*  area_
            ind = wsample(1:it_,w)
            counts[ind,iy_] += 1
            fulldata[k].ix = ind
        end
        if obstype_=="yunknown"  # update ix and iy, consinder i in it_..nbinx, j in 1..nbiny
            w = [θ[iter][i,j] for i in it_:m for j in 1:n] .* area_
            ind = wsample(fulldata[k].indpairs,w)
            counts[ind[1],ind[2]] += 1
            fulldata[k].ix = ind[1]
            fulldata[k].iy = ind[2]
        end
    end
    if mod(iter,50)==0
        println(iter)
    end
end

# compute average of weights
θpostmean_dir = mat2vec(mean(θ))

plotting(θpostmean_dir,truedatagen,binarea,binx,biny,"Dirichlet","contourDirichlet")

traceplotswanted = false
if traceplotswanted
    θiterates = [Any[θ[iter][i,j], iter, "[$i,$j]"] for iter in 1:ITER, i in 1:m, j in 1:n][:]
    iterates_df = DataFrame(w=extractind(θiterates,1), iterate=extractind(θiterates,2), binID=extractind(θiterates,3))
    @rput iterates_df
    R"""
    library(ggplot2)
    library(tidyverse)
    iterates_df %>%  ggplot() + geom_line(aes(x=iterate,y=w)) + facet_wrap(~binID,scales="free") #+theme_minimal()
    """
end
