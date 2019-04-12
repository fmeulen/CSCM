
"""
Evaluate error in binning by piecewise constant density assumption.
"""
function errorftrue_binning(x,y,truedatagen,binx,biny)
    i = indbin(x,binx)
    j = indbin(y,biny)
    area = (biny[j+1]-biny[j]) * (binx[i+1] - binx[i])
    bindens = (Ftrue(binx[i+1],biny[j+1], truedatagen) -
                Ftrue(binx[i+1],biny[j], truedatagen) -
                Ftrue(binx[i],biny[j+1], truedatagen) +
                Ftrue(binx[i],biny[j], truedatagen) )/area
    ftrue(x,y, truedatagen) - bindens
end


"""
Evaluate cdf of piecewise constant density in point (x,y), when the values on the bins are contained in the matrix w, and
the grid specifying the bins is given by binx (horizontal direction) and biny (vertical direction)

Implementation can be made more efficient, but probably not worth the effort.
"""
function cdf_(x,y,w,binx,biny)  # weights correspond to density in bins
    ix = indbin(x,binx)
    iy = indbin(y,biny)
    out = (x-binx[ix])*(y-biny[iy])*w[ix,iy]
    for i in 1:ix-1
        out += (binx[i+1]-binx[i]) * (y-biny[iy])*w[i,iy]
    end
    for j in 1:iy-1
        out += (x-binx[ix]) * (biny[j+1]-biny[j])*w[ix,j]
    end
    for i in 1:ix-1, j in 1:iy-1
        out += (binx[i+1]-binx[i])*(biny[j+1]-biny[j])*w[i,j]
    end
    out
end



function contourplot(weights, truedatagen; titel="Contour plot", gridN=200,binx=binx,biny=biny) # weights correspond to density in bins
    minx, maxx = extrema(binx)
    miny, maxy = extrema(biny)
    gridx = range(minx,stop=maxx-0.001,length=gridN)
    gridy = range(miny,stop=maxy-0.001,length=gridN)
     # gridF =[cdf_(gridx[i],gridy[j],weights,binx,biny) for j in 1:gridN for i in 1:gridN]
     # gridFᵒ = [Ftrue(gridx[i],gridy[j],truedatagen)  for j in 1:gridN for i in 1:gridN]

    gridF =[cdf_(gridx[i],gridy[j],weights,binx,biny)  for i in 1:gridN  for j in 1:gridN]
    gridFᵒ = [Ftrue(gridx[i],gridy[j],truedatagen)   for i in 1:gridN  for j in 1:gridN] # this works fine

    xxgrid = repeat(gridx,inner=gridN)
    yygrid = repeat(gridy,outer=gridN)
    leveldf = DataFrame(x=xxgrid,y=yygrid, posteriormean=gridF, datagenerating=gridFᵒ)

    @rput leveldf
    @rput titel
    R"""
        library(ggplot2)
        library(tidyverse)
        b<- c(0.01,0.05,seq(0.1:1.0,by=0.1))

        d <- leveldf %>% gather(key= curveID, value=z, posteriormean, datagenerating)
        d %>% ggplot(aes(x=x,y=y,z=z,colour=curveID,linetype=curveID)) + stat_contour(breaks=b, size=1.5) +
        theme(legend.position="bottom")+ggtitle(titel) #+theme_minimal()
        ggsave(paste0(titel,".pdf"))
    """
end


function plotapproxerror(truedatagen,binarea,binx,biny;gridN=200)
    # Asses error by binning for true pdf
    minx, maxx = extrema(binx)
    miny, maxy = extrema(biny)
    gridx = range(minx,stop=maxx-0.001,length=gridN)
    gridy = range(miny,stop=maxy-0.001,length=gridN)
    err = Float64[]
    for i in eachindex(gridx)
        for j in eachindex(gridy)
            push!(err, errorftrue_binning(gridx[i],gridy[j],truedatagen,binx,biny))
        end
    end
    d2 = DataFrame(error=err, x =repeat(gridx,inner=gridN), y=repeat(gridy,outer=gridN))
    @rput d2
    R"""
        library(ggplot2)
        library(tidyverse)
        ggplot(data=d2,aes(x, y, fill=error)) + geom_tile() + ggtitle("approximation error")
        ggsave("approxerror.pdf")
    """
end
