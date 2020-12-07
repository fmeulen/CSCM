struct obs{S<:Real, T<:Int} #FIXME need not be mutable.
    t::S
    ix::T
    iy::T
    it::T
    obstype::String
    area::Vector{Float64}  # keep track of bin areas or lengths (in case y is known), needed for updating latent data
    indpairs::Array{Array{T,1},1}  # only needed for case y is unknown
end



function sample_dir(t,ind_yknown, y, (binx, biny), ITER; priorscale = 0.1)
    nsample = length(t)
    m = length(binx) - 1
    n = length(biny) - 1
    counts = zeros(Int64,m,n)  # adjust at each iteration
    fulldata = Vector{obs}(undef,nsample )  # adjust at each iteration

    for k ∈ eachindex(t)
        it = indbin(t[k],binx)
        if k in ind_yknown
            ix = rand(1:it)
            iy = indbin(y[k],biny)
            obstype="yknown"
            area = [ min(binx[i+1],t[k])-binx[i] for i in 1:it]
            indpairs = [ [0,0]] #0
        else
            ix = rand(it:m)
            iy = rand(1:n)
            obstype="yunknown"
            area = [(binx[i+1]-max(t[k],binx[i])) * (biny[j+1] - biny[j])  for i in it:m for j in 1:n]
            indpairs =[ [i,j] for i in it:m for j in 1:n]
        end
        counts[ix,iy] +=1
        fulldata[k] = obs(t[k],ix,iy,it,obstype,area,indpairs)
    end

    priorθ = priorscale *  ones(m,n)

    θ = Vector{Matrix{Float64}}(undef,ITER) # save in each iteration

    ####### Data augmentation algorithm
    for iter ∈ 1:ITER
        # Update weights
        θ[iter] = vec2mat( rand(Dirichlet(vec(counts+priorθ))) , m, n)

        # Update latent data
        for k ∈ eachindex(t)
            ix_ = fulldata[k].ix;  iy_ = fulldata[k].iy; it_ = fulldata[k].it;
            t_ = fulldata[k].t
            obstype_ = fulldata[k].obstype
            area_ = fulldata[k].area

            counts[ix_,iy_] += -1
            if obstype_=="yknown"  # update ix, consider i in 1..it_,  j=iy_
                w = [θ[iter][i,iy_] for i in 1:it_] .*  area_
                ind = wsample(1:it_,w)
                counts[ind,iy_] += 1
                fulldata = @set fulldata[k].ix = ind
            end
            if obstype_=="yunknown"  # update ix and iy, consinder i in it_..nbinx, j in 1..nbiny
                w = [θ[iter][i,j] for i in it_:m for j in 1:n] .* area_
                ind = wsample(fulldata[k].indpairs,w)
                counts[ind[1],ind[2]] += 1
                fulldata = @set fulldata[k].ix = ind[1]
                fulldata = @set fulldata[k].iy = ind[2]
            end
        end
        if mod(iter,250)==0
            println(iter)
        end
    end

    # compute average of weights
    θ
end
