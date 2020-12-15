# code to sample from the posterior for both the Dirichlet and logistic-Normal-graphLaplacina prior

graphlaplacian(m,n; pow=1) = (Matrix(lap(grid2(m,n))) + I/(m*n)^2)^pow

"""
    dirichlet(ci, bins::Bins, IT, Πτ; τinit = 1.0, δ=0.1, printskip=5000)

ci:: CensoringInfo (contains info on the observations)
bins:: Bins (contains info on the bins)
IT: number of iterations
Πτ: prior distribution on τ

Returns:
θsave, τsave, acc
"""
function dirichlet(ci, bins::Bins, IT, Πτ; τinit = 1.0, δ=0.1, printskip=5000)
    m, n = bins.m, bins.n
    N = m*n
    # initialise counts of bins by random assignment
    counts_fulldata = zeros(N)
    for i ∈ eachindex(ci)
        loc = rand(ci[i].ind)
        counts_fulldata[loc] += 1
    end

    τ = τinit
    θ = rand(Dirichlet(counts_fulldata .+ τ))

    θsave = zeros(IT,N)
    τsave = zeros(IT)
    τsave[1] = τ
    θsave[1,:] = θ
    acc = 0

    for it ∈ 2:IT
        # update counts
        counts_fulldata .= zeros(N)
        for i ∈ eachindex(ci)
            γ = ci[i]
            loc = wsample(γ.ind, θ[γ.ind] .* γ.fracarea)
            counts_fulldata[loc] += 1
        end
        # update θ
        rand!(Dirichlet( counts_fulldata .+ τ), θ)
        θsave[it,:]  = θ
        # update τ
        τᵒ = τ * exp(δ*randn())
        A = logpdf(Dirichlet(N,τᵒ), θ)  -
            logpdf(Dirichlet(N,τ), θ) +
            logpdf(Πτ, τᵒ) - logpdf(Πτ, τ) +
            log(τᵒ) - log(τ)
        if log(rand()) < A
            τ = τᵒ
            acc += 1
        end
        τsave[it] = τ

        if mod(it,printskip)==0
            frac_accepted = round.(acc/it;digits=2)
            @show (it, frac_accepted)
        end
    end
    θsave, τsave, acc
end



"""
    loglik!(θ, τ, z,  Uinv, ci, Πτ)

Computes θ and the sum of the logprior and loglikelihood for (τ,z)
θ is written into

Uinv = inv(L.chol.U), where L is the graphLaplacian matrix
ci::CensoringInfo contains information on the observations
"""
function loglik!(θ, τ, z,  Uinv, ci, Πτ)
    softmax!(θ, Uinv * z * √τ)
    ll = logpdf(Πτ, τ)
    @inbounds for i ∈ eachindex(ci)
        γ = ci[i]
        ll += log(dot(θ[γ.ind], γ.fracarea))
    end
    θ, ll
end

"""
    pcn(ci, bins::Bins, IT, Πτ; ρ = 0.95, τinit = 1.0, δ=0.1, printskip=5000)

ci:: CensoringInfo (contains info on the observations)
bins:: Bins (contains info on the bins)
IT: number of iterations
Πτ: prior distribution on τ

Returns:
θsave, τsave, acc, ρ
"""
function pcn(ci, bins::Bins, IT, Πτ; ρ = 0.95, τinit = 1.0, δ=0.1, printskip=5000)
    m, n = bins.m, bins.n
    L = PDMat(graphlaplacian(m,n))
    Uinv = inv(L.chol.U)

    N = m*n
    z = randn(N)
    τ = τinit
    θ, ll = loglik!(zeros(N), τ, z, Uinv, ci, Πτ)

    # cache arrays
    zᵒ = zeros(N)
    znew = zeros(N)
    θᵒ = zeros(N)

    ρc = sqrt(1.0 - ρ^2)
    θsave = zeros(IT,N)
    τsave = zeros(IT)
    θsave[1,:]  = θ
    τsave[1] = τ
    acc = [1, 1]
    for it ∈ 2:IT
        # update z
        randn!(znew)
        zᵒ .= ρ * z + ρc * znew
        θᵒ, llᵒ = loglik!(θᵒ, τ, zᵒ, Uinv, ci, Πτ)
        if log(rand()) < llᵒ-ll
            z .= zᵒ
            θ .= θᵒ
            ll = llᵒ
            acc[1] += 1
        end
        θsave[it,:]  = θ

        # update τ
        τᵒ = τ * exp(δ * randn())
        θᵒ, llᵒ = loglik!(θᵒ, τᵒ, z, Uinv, ci, Πτ)
        A = llᵒ - ll + log(τᵒ) - log(τ)
        if log(rand()) < A
            τ = τᵒ
            θ .= θᵒ
            ll = llᵒ
            acc[2] += 1
        end
        τsave[it] = τ

        if mod(it,printskip)==0
            frac_acc = round.(acc/it;digits=2)
            @show (it, frac_acc)
        end
    end
    @show "Fraction of accepted pCN steps equals: " acc/IT
    θsave, τsave, acc, ρ
end
