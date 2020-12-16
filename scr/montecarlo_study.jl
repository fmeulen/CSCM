
dist = Mixture(0.7,X2plusyRev(), X2plusy());
m = 25; n= 50
#dist = Xplusy(); m = 25; n= 25
nsample = [100, 250, 500]

Πdir = Exponential(1.0)   #Uniform(0.01, 1.0)
Πτ = Exponential(1.0)
IT = 20_000 # nr of iterations for Dirichlet prior
nMC = 100

function mcmcstudy(dist, nsample::Vector, m, n, Πdir, Πτ, nMC, IT)
    BI = div(IT,3) # nr of burnin iters
    bi = BI:IT
    bins = Bins(dist, m, n)
    θ0, xx, yy = binprob(dist,bins)

    outdir = Array{Float64,1}[]
    outgl = Array{Float64,1}[]
    for k ∈ eachindex(nsample)
        𝒲dir = zeros(nMC)
        𝒲gl = zeros(nMC)
        for j ∈ 1:nMC
            x, y, t, ind_yknown, ind_yunknown = gencensdata(dist, nsample[k])
            ci = construct_censoringinfo(t, y, ind_yknown, ind_yunknown, bins)
            θdir, τdir, accdir, iters_saved = dirichlet(ci, bins, IT, Πdirτ)
            θgl, τgl, accgl, ρ, iters_saved = pcn(ci, bins, IT, Πτ; ρ=.96, δ=0.6)
            θ̄dir = [mean(x[bi]) for x ∈ eachcol(θdir)]
            θ̄gl = [mean(x[bi]) for x ∈ eachcol(θgl)]
            𝒲dir_,  𝒲gl_ = wasserstein(θ̄dir, θ̄gl, θ0, bins)
            𝒲dir[j] = 𝒲dir_
            𝒲gl[j] = 𝒲gl_
        end
        push!(outdir, 𝒲dir)
        push!(outgl, 𝒲gl)
    end
    outdir, outgl
end

outdir, outgl = mcmcstudy(dist, nsample, m, n, Πdir, Πτ, nMC, 100)

out = DataFrame(D = vec(hcat(outdir...)), LNGL = vec(hcat(outgl...)), samplesize=repeat(nsample,inner = nMC))
odir = mkpath("./out/mcstudy")
CSV.write(joinpath(odir,"mcstudy.csv"),out)
