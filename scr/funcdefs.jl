"""
Conversion vector and matrix format
Example usage
    a =  reshape(1:6,2,3)
    vec2mat(mat2vec(a),2,3)-a
"""
mat2vec(x) = vec(x) # colunmwise filling
vec2mat(x,m,n) = reshape(x,m,n)

extractind(w,i) = map(x -> x[i], w)

# look up in which bin y is located using
indbin(y, bin) = (y >= bin[end]) ? error("bin does not exist") : findfirst(x -> x > y, bin) - 1

abstract type GeneratingDistribution end

struct Xplusy <: GeneratingDistribution  end
struct XplusyRev <: GeneratingDistribution end
struct Uniform2D <: GeneratingDistribution end
struct X2plusy <: GeneratingDistribution end
struct X2plusyRev <: GeneratingDistribution end
struct Mixture <: GeneratingDistribution
    p::Float64
    comp1::GeneratingDistribution
    comp2::GeneratingDistribution
end

loss(x, c) = (x-c)^2
err(dist::GeneratingDistribution, c) =  x -> loss(dens(dist,x), c)

#-------------
function gendata(::Xplusy, N)
        x = -0.5 .+ 0.5 *sqrt.(1 .+ 8*rand(N))
        y = -x .+ sqrt.(x.^2 .+ (2x .+ 1) .* rand(N))
        x, y
end
binprob(::Xplusy, xm, x, ym, y) =  0.5*(y-ym)*(x^2-xm^2) + 0.5*(x-xm)*(y^2-ym^2)
support(::Xplusy) = [(0.0, 1.0), (0.0, 1.0)]
dens(::Xplusy, x) = x[1] + x[2]
#-------------
function gendata(dist::XplusyRev,N)
    x, y = gendata(Xplusy(), N)
    1.0 .- x, 1.0 .- y
end
binprob(::XplusyRev, xm, x, ym, y) = binprob(Xplusy(), 1.0 - x, 1.0 - xm, 1.0 - y, 1.0 - ym)
support(::XplusyRev) = support(Xplusy())
dens(::XplusyRev,x) = 100.0 #@error "not implemented yet"

#-------------
function gendata(::Uniform2D, N)
    rand(N), rand(N)
end
binprob(::Uniform2D, xm, x, ym, y) = (x-xm)*(y-ym)
support(::Uniform2D) = [(0.0, 1.0), (0.0, 1.0)]
dens(::Uniform2D, x) = 1.0
#-------------
function gendata(::X2plusy, N)
    u = rand(N)
    x = (2u .+ sqrt.(4*u.^2 .+ 1)).^(1/3) .- abs.(2u .- sqrt.(4*u.^2 .+ 1)).^(1/3)
    y = -x.^2 + sqrt.(x.^4 .+ 4*(x.^2 .+ 1) .* rand(N))
    x, y
end
binprob(::X2plusy, xm, x, ym, y) = (3/8)*((y-ym)*(x^3 - xm^3)/3 + 0.5*(x-xm)*(y^2-ym^2))
support(::X2plusy) = [(0.0, 1.0), (0.0, 2.0)]
dens(::X2plusy, x) = (3/8)*(x[1]^2 + x[2])
#-------------
function gendata(dist::X2plusyRev,N)
    x, y = gendata(X2plusy(), N)
    1.0 .- x,  y
end
binprob(::X2plusyRev, xm, x, ym, y) = binprob(X2plusy(), 1.0 - x, 1.0 - xm, ym, y)
support(::X2plusyRev) = support(X2plusy())
dens(::X2plusyRev,x) = 100.0 #@error "not implemented yet"
#-------------
function gendata(dist::Mixture, N)
    u = rand(Bernoulli(dist.p), N)
    x1, y1 = gendata(dist.comp1,N)
    x2, y2 = gendata(dist.comp2, N)
    u .* x1 + (1 .- u) .* x2, u .* y1 + (1 .- u) .*y2
end
binprob(dist::Mixture, xm, x, ym, y) = dist.p * binprob(dist.comp1, xm, x, ym, y) + (1.0 - dist.p) * binprob(dist.comp2, xm, x, ym, y)
support(dist::Mixture) = support(dist.comp1)  ## FIXME
dens(dist::Mixture, x) = dist.p * dens(dist.comp1, x) +(1.0 - dist.p) * dens(dist.comp2, x)




function gencensdata(dist::GeneratingDistribution, N)
    x, y = gendata(dist,N)
    t = sqrt.(rand(length(x)))
    # find indices of case where x<t (y observed) and x>= t (y not observed)
    ind_yknown = findall(x.<t)
    ind_yunknown = findall(x.>=t)
    # So as available data we are given ind_yknown, ind_yunknown, t and y[ind_yknown]
    x, y, t, ind_yknown, ind_yunknown
end


# construct bins
struct Bins{T}
    m::Int64
    n::Int64
    binx::T
    biny::T
end

function Bins(dist::GeneratingDistribution, m, n)
    𝒮 = support(dist)
    binx = range(𝒮[1][1], stop=𝒮[1][2], length=m+1)
    biny = range(𝒮[2][1], stop=𝒮[2][2], length=n+1)
    Bins(m, n, binx, biny)
end




"""
    binprob(dist::GeneratingDistribution,bins)

returns xvals, yvals and binvals, where binsvals is computed for distribution `dist`
"""
function binprob(dist::GeneratingDistribution,bins)
    m, n, binx, biny = bins.m, bins.n, bins.binx, bins.biny
    xx = repeat(binx[2:end],inner=n)
    yy = repeat(biny[2:end],outer=m)
    out = Float64[]
    for i in 1:m,  j in 1:n
        push!(out, binprob(dist, binx[i], binx[i+1], biny[j], biny[j+1]))
    end
    out, xx, yy
end

function binerror(dist::GeneratingDistribution, bins, θ)
    m, n, binx, biny = bins.m, bins.n, bins.binx, bins.biny
    c = vec2mat(θ,m,n)
    out = Float64[]
    for i in 1:m,  j in 1:n
        area = (binx[i+1]-binx[i]) * (biny[j+1] - biny[j])
        push!(out, hcubature(err(dist,c[i,j]/area), [binx[i], biny[j]], [binx[i+1], biny[j+1]])[1] )
    end
    out
end

struct CensoringInfo{S<:Number, T<:Number}
    fracarea::Vector{S}         # keep track of fraction of bin areas
    ind::Vector{T}              # corresponding indices
end

function construct_censoringinfo(t, y, ind_yknown, ind_yunknown, bins)
    m, n, binx, biny = bins.m, bins.n, bins.binx, bins.biny
    nsample = length(t)
    # construct censoringinfo
    ci = Vector{CensoringInfo}(undef,nsample)
    for k ∈ ind_yknown
        it = indbin(t[k],binx)
        iy = indbin(y[k],biny)
        fa =  [ (min(binx[i+1],t[k])-binx[i])/(binx[i+1]-binx[i]) for i ∈ 1:it]
        ind = [iy + ℓ*n for ℓ ∈ 0:(it-1)]
        ci[k] = CensoringInfo(fa, ind)
    end
    for k ∈ ind_yunknown
        it = indbin(t[k],binx)
        fa = [(binx[i+1]-max(t[k],binx[i]))/(binx[i+1] - binx[i])  for i ∈ it:m for j ∈ 1:n]
        ind = collect(((it-1)*n+1):(m*n))
        ci[k] = CensoringInfo(fa,ind)
    end
    ci
end




function wasserstein(θ̄dir, θ̄gl, θ0, bins::Bins; p=1)
    points_x = [mean(bins.binx[i:i+1]) for i ∈ 1:bins.m]
    points_y = [mean(bins.biny[i:i+1]) for i ∈ 1:bins.n]
    out = [[u, v] for u in points_x for v in points_y]
    coordmat = zeros(bins.m*bins.n, 2)
    for i in eachindex(out)
        coordmat[i,:] = out[i]
    end

    thgl = θ̄gl
    thdir = θ̄dir
    th0 = θ0
    @rput thgl thdir th0 m n coordmat p
    R"""
    library(transport)
    gl = wpp(coordmat, thgl)
    dir = wpp(coordmat, thdir)
    tr0 = wpp(coordmat, th0)
    was_dir = wasserstein(dir, tr0, p=p)
    was_gl = wasserstein(gl, tr0, p=p)
    """
    @rget was_gl was_dir
    was_dir, was_gl
end
