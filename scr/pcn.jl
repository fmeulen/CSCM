


function loglik!(θ, τ, z,  L, ci)
    softmax!(θ, √(τ) * (L*z))
    ll = 0.0
    @inbounds for i ∈ eachindex(ci)
        γ = ci[i]
        ll += log(dot(θ[γ.ind], γ.fracarea))
    end
    θ, ll
end

function pcn(ci, IT; ρ=0.9, τ = 1.0)
    L = graphlaplacian(m,n)
    Lchol = PDMat(L).chol.L

    N = m*n
    z = randn(N)

    θ, ll = loglik!(zeros(N), τ, z, Lchol, ci)

    # cache arrays
    zᵒ = zeros(N)
    znew = zeros(N)
    θᵒ = zeros(N)


    ρc = sqrt(1.0-ρ^2)
    θsave = zeros(IT,N)
    θsave[1,:]  = θ
    acc = 0
    for i ∈ 2:IT
        randn!(znew)
        zᵒ .= ρ * z + ρc * znew
        θᵒ, llᵒ = loglik!(θᵒ, τ, zᵒ, Lchol, ci)
        if log(rand()) < llᵒ-ll
            z .= zᵒ
            θ .= θᵒ
            ll = llᵒ
            acc += 1
        end
        #push!(θsave, θ)
        θsave[i,:]  = θ
    end
    θsave, acc
end

IT = 10_000; BI = div(IT,2)
ci = construct_censoringinfo(t, (binx,biny), ind_yknown, ind_yunknown)
θsave, acc = pcn(ci,IT; τ=0.1)
θ̄gl = [mean(x[BI:end]) for x ∈ eachcol(θsave)]
acc/IT

Plots.plot(θsave[:,3])
Plots.plot(θsave[:,50])
