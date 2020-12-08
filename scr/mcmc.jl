# code to sample from the posterior for both the Dirichlet and logistic-Normal-graphLaplacina prior

struct CensoringInfo{S<:Number, T<:Number}
    fracarea::Vector{S}         # keep track of fraction of bin areas
    ind::Vector{T}              # corresponding indices
end

function construct_censoringinfo(t, y, (binx,biny), ind_yknown, ind_yunknown)
    nsample = length(t)
    m = length(binx) - 1
    n = length(biny) - 1
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

graphlaplacian(m,n) = Matrix(lap(grid2(m,n))) + I/(m*n)^2

"""
    dirichlet(ci, (m, n), IT; τinit = 1.0, δ=0.1, priorτ = InverseGamma(0.1,0.1))

ci:: CensoringInfo
(m, n): number of horizontal and vertical bins
IT: number of iterations

Returns
θsave, τsave, acc
"""

function dirichlet(ci, (m, n), IT;
            τinit = 1.0, δ=0.1, priorτ = InverseGamma(0.1,0.1))
    N = m*n
    # initialise counts of bins by random assignment
    counts_fulldata = zeros(N)
    for i ∈ eachindex(ci)
        loc = rand(ci[i].ind)
        counts_fulldata[loc] += 1
    end
    τ = τinit
    θsave = zeros(IT,N)
    τsave = zeros(IT)
    θ = rand(Dirichlet( counts_fulldata .+ τ))
    τsave[1] = τ
    θsave[1,:]  = θ
    acc = 0
    for it ∈ 2:IT
        # update counts
        for i ∈ eachindex(ci)
            γ = ci[i]
            loc = wsample(γ.ind, θ[γ.ind] .* γ.fracarea)
            counts_fulldata[loc] += 1
        end
        # update θ
        θ = rand(Dirichlet( counts_fulldata .+ τ))
        θsave[it,:]  = θ
        # update τ
        τᵒ = τ * exp(δ*randn())
        A = logpdf(Dirichlet(N,τᵒ), θ)  -
            logpdf(Dirichlet(N,τ), θ) +
            logpdf(priorτ, τᵒ) - logpdf(priorτ, τ) +
            log(τᵒ) - log(τ)
        if log(rand()) < A
            τ = τᵒ
            acc += 1
        end
        τsave[it] = τ


    end
    θsave, τsave, acc
end



"""
    loglik!(θ, τ, z,  Uinv, ci)

Computes θ and loglikelihood for (τ,z)
θ is written into

Uinv = inv(L.chol.U), where L is the graphLaplacian matrix
ci::CensoringInfo contains information on the observations
"""
function loglik!(θ, τ, z,  Uinv, ci)
    softmax!(θ, Uinv * z * √τ)
    ll = 0.0
    @inbounds for i ∈ eachindex(ci)
        γ = ci[i]
        ll += log(dot(θ[γ.ind], γ.fracarea))
    end
    θ, ll
end

"""
    pcn(ci, (m, n), IT; ρ = 0.95, τinit = 1.0, δ=0.1, priorτ = InverseGamma(0.1,0.1))

ci:: CensoringInfo
(m, n): number of horizontal and vertical bins
IT: number of iterations

Returns
θsave, τsave, acc, ρ
"""
function pcn(ci, (m, n), IT; ρ = 0.95, τinit = 1.0, δ=0.1, priorτ = InverseGamma(0.1,0.1))

    L = PDMat(graphlaplacian(m,n))
    Uinv = inv(L.chol.U)

    N = m*n
    z = randn(N)
    τ = τinit
    θ, ll = loglik!(zeros(N), τ, z, Uinv, ci)

    # cache arrays
    zᵒ = zeros(N)
    znew = zeros(N)
    θᵒ = zeros(N)

    ρc = sqrt(1.0-ρ^2)
    θsave = zeros(IT,N)
    τsave = zeros(IT)
    θsave[1,:]  = θ
    τsave[1] = τ
    acc = [1, 1]
    for i ∈ 2:IT
        # update z
        randn!(znew)
        zᵒ .= ρ * z + ρc * znew
        θᵒ, llᵒ = loglik!(θᵒ, τ, zᵒ, Uinv, ci)
        if log(rand()) < llᵒ-ll
            z .= zᵒ
            θ .= θᵒ
            ll = llᵒ
            acc[1] += 1
        end
        θsave[i,:]  = θ

        # update τ
        τᵒ = τ * exp(δ*randn())
        θᵒ, llᵒ = loglik!(θᵒ, τᵒ, z, Uinv, ci)
        A = sum(llᵒ) - sum(ll) +
            logpdf(priorτ, τᵒ) - logpdf(priorτ, τ) +
            log(τᵒ) - log(τ)
        if log(rand()) < A
            τ = τᵒ
            θ .= θᵒ
            ll = llᵒ
            acc[2] += 1
        end
        τsave[i] = τ

        if mod(i,1000)==0
            frac_acc = round.(acc/i;digits=2)
            @show (i, frac_acc)
        end

    end
    @show "Fraction of accepted pCN steps equals: $(acc/IT)"
    θsave, τsave, acc, ρ
end
