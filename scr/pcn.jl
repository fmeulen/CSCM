function loglik!(θ, τ, z,  U, ci)
    softmax!(θ, U * z/√τ)
    ll = 0.0
    @inbounds for i ∈ eachindex(ci)
        γ = ci[i]
        ll += log(dot(θ[γ.ind], γ.fracarea))
    end
    θ, ll
end

function pcn(t,ind_yknown, y, (binx, biny), IT; ρ = 0.95, τinit = 1.0)
    ci = construct_censoringinfo(t, (binx,biny), ind_yknown, ind_yunknown)
    m, n = length(binx) - 1, length(biny) - 1
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
    acc = 0
    for i ∈ 2:IT
        randn!(znew)
        zᵒ .= ρ * z + ρc * znew
        θᵒ, llᵒ = loglik!(θᵒ, τ, zᵒ, Uinv, ci)
        if log(rand()) < llᵒ-ll
            z .= zᵒ
            θ .= θᵒ
            ll = llᵒ
            acc += 1
        end
        θsave[i,:]  = θ
        τ = rand(InverseGamma(0.5N + 0.1, 0.5PDMats.quad(L, z) + 0.1))
        τsave[i] = τ
    end
    @show "Fraction of accepted pCN steps equals: $(acc/IT)"
    θsave, τsave, acc, ρ
end
