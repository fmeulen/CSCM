"""
Conversion vector and matrix format
Example usage
    a =  reshape(1:6,2,3)
    vec2mat(mat2vec(a),2,3)-a
"""
mat2vec(x) = vec(x) # colunmwise filling
vec2mat(x,m,n) = reshape(x,m,n)

extractind(w,i) = map(x -> x[i],w)

indbin(x_,binx) = (x_ >= binx[end]) ? error("bin does not exist") : findfirst(x -> x > x_, binx) - 1  # look up in which bin x_ is located



abstract type GeneratingDistribution end

struct Xplusy <: GeneratingDistribution  end
struct XplusyRev <: GeneratingDistribution end
struct Uniform2D <: GeneratingDistribution end
struct X2plusy <: GeneratingDistribution end
struct Mixture <: GeneratingDistribution
    p::Float64
    comp1::GeneratingDistribution
    comp2::GeneratingDistribution
end

loss(x, c) = (x-c)^2
err(dist::GeneratingDistribution, c) =  x -> loss(dens(dist,x), c)

function gendata(::Xplusy, N)
        x = -0.5 .+ 0.5 *sqrt.(1 .+ 8*rand(N))
        y = -x .+ sqrt.(x.^2 .+ (2x .+ 1) .* rand(N))
        x, y
end
binprob(::Xplusy, xm, x, ym, y) =  0.5*(y-ym)*(x^2-xm^2) + 0.5*(x-xm)*(y^2-ym^2)
support(::Xplusy) = [(0.0, 1.0), (0.0, 1.0)]
dens(::Xplusy, x) = x[1] + x[2]


function gendata(::Uniform2D, N)
    rand(N), rand(N)
end
binprob(::Uniform2D, xm, x, ym, y) = (x-xm)*(y-ym)
support(::Uniform2D) = [(0.0, 1.0), (0.0, 1.0)]
dens(::Uniform2D, x) = 1.0

function gendata(::X2plusy, N)
    u = rand(N)
    x = (2u .+ sqrt.(4*u.^2 .+ 1)).^(1/3) .- abs.(2u .- sqrt.(4*u.^2 .+ 1)).^(1/3)
    y = -x.^2 + sqrt.(x.^4 .+ 4*(x.^2 .+ 1) .* rand(N))
    x, y
end
binprob(::X2plusy, xm, x, ym, y) = (3/8)*((y-ym)*(x^3 - xm^3)/3 + 0.5*(x-xm)*(y^2-ym^2))
support(::X2plusy) = [(0.0, 1.0), (0.0, 2.0)]
dens(::X2plusy, x) = (3/8)*(x[1]^2 + x[2])


function gendata(dist::Mixture, N)
    u = rand(Bernoulli(dist.p), N)
    x1, y1 = gendata(dist.comp1,N)
    x2, y2 = gendata(dist.comp2, N)
    u .* x1 + (1 .- u) .* x2, u .* y1 + (1 .- u) .*y2
end
binprob(dist::Mixture, xm, x, ym, y) = dist.p * binprob(dist.comp1, xm, x, ym, y) + (1.0 - dist.p) * binprob(dist.comp2, xm, x, ym, y)
support(::Mixture) = [(0.0, 1.0), (0.0, 1.0)] #FIXME
dens(::Mixture, x) = (3/8)*(x[1]^2 + x[2]) # FIXME

function gendata(dist::XplusyRev,N)
    x, y = gendata(Xplusy(), N)
    1.0 .- x, 1.0 .- y
end
binprob(::XplusyRev, xm, x, ym, y) = binprob(Xplusy(), 1.0 - x, 1.0 - xm, 1.0 - y, 1.0 - ym)
support(::XplusyRev) = [(0.0, 1.0), (0.0, 1.0)]
dens(::XplusyRev) = @error "not implemented yet"

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



#
# """
#     gendata(truedatagen, N)
#
# Generate data for current status continuous mark model.
# """
# function gendata(truedatagen, N; Ntest=10)
#     t = sqrt.(rand(N))
#     if truedatagen=="x+y"
#         x = -0.5 .+ 0.5 *sqrt.(1 .+ 8*rand(N))
#         y = -x .+ sqrt.(x.^2 .+ (2x .+ 1) .* rand(N))
#     end
#     if truedatagen=="uniform"
#         x = rand(N)
#         y = rand(N)
#     end
#     if truedatagen=="(3/8)(x2+y)"
#         u = rand(N)
#         x = (2u .+ sqrt.(4*u.^2 .+ 1)).^(1/3) .- abs.(2u .- sqrt.(4*u.^2 .+ 1)).^(1/3)
#         y = -x.^2 + sqrt.(x.^4 .+ 4*(x.^2 .+ 1) .* rand(N))
#     end
#
#
#     if truedatagen=="testcase"
#         bx_ = range(0.0, stop=1.0, length=Ntest)
#         by_ = range(0.0, stop=1.0, length=Ntest)
#         θ0 = zeros(Ntest,Ntest)
#         N1 = MvNormal([1., 2.]/ 3.0, 0.1*PDMat([0.2 0.0 ; 0.0 0.2]))
#         N2 = MvNormal([3., 1.]/ 4.0, 0.1*PDMat([0.2 0.0 ; 0.0 0.2]))
#         for i in 1:Ntest
#             for j in 1:Ntest
#                 u = [bx_[i], by_[j]]
#                 θ0[i,j] =  0.4 * pdf(N1,u) + 0.6 * pdf(N2,u)
#             end
#         end
#         w = wsample(1:(Ntest^2), vec(θ0) , N, replace=true)
#         vals = collect(1:Ntest)/Ntest
#         x = zeros(N)
#         y = zeros(N)
#         for i in 1:N
#             @show w[i]
#             ind = div(w[i],Ntest)
#             x[i] = vals[1+ind] + 0.01* randn()
#             y[i] = vals[1 + mod(w[i],Ntest)] + 0.01* randn()
#         end
#     end
#     # find indices of case where x<t (y observed) and x>= t (y not observed)
#     ind_yknown = findall(x.<t)
#     ind_yunknown = findall(x.>=t)
#     # So as available data we are given ind_yknown, ind_yunknown, t and y[ind_yknown]
#
#     x, y, t, ind_yknown, ind_yunknown
# end
#
# function Ftrue(x,y,truedatagen)
#     if truedatagen=="uniform"
#         out = x*y
#     end
#     if truedatagen=="x+y"
#         out = 0.5*x^2*y + 0.5*x*y^2
#     end
#     if truedatagen=="(3/8)(x2+y)"
#         out = (3/8)*(y*(x^3)/3 + 0.5*x*y^2)
#     end
#     out
# end

# function ftrue(x,y,truedatagen)
#     if truedatagen=="uniform"
#         out = 1
#     end
#     if truedatagen=="x+y"
#         out = x+y
#     end
#     if truedatagen=="(3/8)(x2+y)"
#         out = (3/8)*(x^2 + y)
#     end
#     out
# end

"""
Compute the true bin probabilities for data generated by 'truedatagen'.
Probabilities are stacked into a vector where the inner loop goes over the y-direction.
"""
# function binprobtrue(dist::GeneratingDistribution,  binx, biny)
#     m, n = length(binx)-1, length(biny)-1
#     out = Float64[]
#     for i in 1:m
#         for j in 1:n
#             val = Ftrue(dist, binx[i+1],biny[j+1]) -
#                         Ftrue(dist, binx[i+1],biny[j]) -
#                         Ftrue(dist, binx[i],biny[j+1]) +
#                         Ftrue(dist, binx[i],biny[j])
#             push!(out, val)
#         end
#     end
#     out
# end
#
#
# """
# Compare bin probabilities in pweights to true bin probabilities computed using true datagen
# """
# function θtrue(dist::GeneratingDistribution,binx,biny)
#     m = length(binx)-1 ; n = length(biny)-1
#     xx = repeat(binx[2:end],inner=n)
#     yy = repeat(biny[2:end],outer=m)
#     θ0 = binprobtrue(dist, binx,biny)
#     θ0, xx, yy
# end




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

function binerror(dist::GeneratingDistribution, bins, θ̄)
    m, n, binx, biny = bins.m, bins.n, bins.binx, bins.biny
    c = vec2mat(θ̄,m,n)
    out = Float64[]
    for i in 1:m,  j in 1:n
            push!(out, hcubature(err(dist,c[i,j]), [binx[i], binx[i+1]], [biny[j], biny[j+1]])[1] )

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
