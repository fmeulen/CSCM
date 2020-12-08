function dirichlet(t, ind_yknown, y, (binx, biny), IT; priorscale = 0.1)
    ci = construct_censoringinfo(t, y, (binx,biny), ind_yknown, ind_yunknown)
    m, n = length(binx) - 1, length(biny) - 1
    N = m*n
    # initialise counts of bins by random assignment
    counts_fulldata = zeros(N)
    for i ∈ eachindex(t)
        loc = rand(ci[i].ind)
        counts_fulldata[loc] += 1
    end
    θsave = zeros(IT,N)
    θ = rand(Dirichlet( counts_fulldata .+ priorscale))
    θsave[1,:]  = θ
    for it ∈ 2:IT
        for i ∈ eachindex(t)
            γ = ci[i]
            loc = wsample(γ.ind, θ[γ.ind] .* γ.fracarea)
            counts_fulldata[loc] += 1
        end
        θ = rand(Dirichlet( counts_fulldata .+ priorscale))
        θsave[it,:]  = θ
    end
    θsave
end
