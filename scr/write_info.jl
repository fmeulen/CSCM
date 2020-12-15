function  processoutput(θdir, τdir, accdir, θgl, τgl, accgl, ρ,
        dist, nsample, bins, BI, IT, outdir)

    bi = BI:IT
    # compute posterior mean
    θ̄dir = [mean(x[bi]) for x ∈ eachcol(θdir)]
    θ̄gl = [mean(x[bi]) for x ∈ eachcol(θgl)]

    # traceplots
    tr1 = Plots.plot(θdir[:,1], label="θ[1] D")
    tr2 = Plots.plot(θdir[:,10], label="θ[10] D")
    tr3 = Plots.plot(log.(τdir), label="log(τ) D")
    tr4 = Plots.plot(θgl[:,1], label="θ[1] LNGL")
    tr5 = Plots.plot(θgl[:,10], label="θ[10] LNGL")
    tr6 = Plots.plot(log.(τgl), label="log(τ) LNGL")
    lay = @layout [a b c; d e f]
    P = Plots.plot(tr1, tr2, tr3, tr4, tr5, tr6, layout=lay)
    savefig(P, joinpath(outdir,"traceplots_jl.pdf"))

    dtrace = DataFrame(iterate=1:IT,theta1=θgl[:,1], theta10=θgl[:,10],
                            theta11=θgl[:,11], logtau=log.(τgl))
    CSV.write("./out/tracepcn.csv",dtrace)

    # true binprobs
    θ0, xx, yy = binprob(dist,bins)
    errdir = binerror(dist, bins, θ̄dir)
    errgl = binerror(dist, bins, θ̄gl)

    labels = repeat(["D", "LNGL", "true"], inner=length(θ0))
    d = DataFrame(binprob = [θ̄dir; θ̄gl; θ0],
        loss = [errdir; errgl; 0.0*errgl], x =[xx; xx; xx], y= [yy; yy; yy],
            method=labels)
    CSV.write(joinpath(outdir,"binprobs.csv"),d)

    distdir = norm(θ0 - θ̄dir,1)
    distgl = norm(θ0 - θ̄gl,1)
    @show distgl/distdir


    # write observations to csv file
    yobserved = fill("yes",nsample)
    yobserved[ind_yunknown] .= "no"
    d = DataFrame(x=x,y=y,t=t,yobserved=yobserved)
    CSV.write(joinpath(outdir,"observations.csv"),d)

    # heatmaps
    xs = [string("x", i) for i = 1:bins.m]
    ys = [string("y", i) for i = 1:bins.n]
    h1 = heatmap(xs, ys, reshape(θ0, bins.n, bins.m),title="true")
    h2 = heatmap(xs, ys, reshape(θ̄dir, bins.n, bins.m),title="D")
    h3 = heatmap(xs, ys, reshape(θ̄gl, bins.n, bins.m),title="LNGL")
    l = @layout [a; b; c]
    PP = plot(h1,h2,h3,layout = l)
    savefig(PP, joinpath(outdir,"heatmaps_jl.pdf"))

    # Wasserstein distance
    𝒲dir,  𝒲gl = wasserstein(θ̄dir, θ̄gl, θ0, bins)
    @printf("𝒲dir = %E \n", 𝒲dir)
    @printf("𝒲gl = %E", 𝒲gl)
    @show  𝒲gl/𝒲dir


    # write info to file

    facc = open(joinpath(outdir,"info.txt"),"w")
        write(facc, "Data choice: ",string(dist),"\n")
        write(facc, "Sample size: ", string(nsample), "\n")
        write(facc, "Fraction in data with mark observed: ",string(length(ind_yknown)/nsample),"\n\n")

        write(facc, "bin info", string(bins), "\n")
        write(facc, "Number of iterations: ",string(IT),"\n")
        write(facc, "Number of burnin iterations: ",string(BI),"\n")

        if false
            write(facc, "Turing Sampler: ", string(sp),"\n")
            write(facc, "Number of iterations for Turing: ",string(ITERgl),"\n")
            write(facc, "Number of burnin iterations for Turing: ",string(BIgl),"\n")
        end

        # write(facc, "Wasserstein distance binprobs Dirichlet: ", string(wdir),"\n")
        # write(facc, "Wasserstein distance binprobs Graph Laplacian: ", string(wgl),"\n")
        # write(facc, "Ratio: ", string(ratio), "\n\n")

        write(facc, "Parameter rho for pCN:", string(ρ), "\n")
        write(facc, "Fraction of accepted pCN steps: ", string(accgl[1]/IT),"\n\n")
        write(facc, "Fraction of accepted τ-update steps for LNGL: ", string(accgl[2]/IT),"\n")

        write(facc, "Fraction of accepted τ-update steps for Dir prior: ", string(accdir/IT),"\n")

        write(facc, "Wasserstein distances: \n")
        write(facc, "for Dir-prior", 𝒲dir, "\n")
        write(facc, "for LNGL-prior", 𝒲gl, "\n")
        write(facc, "ratio:", round(𝒲gl/𝒲dir; digits=3),"\n\n")

    close(facc)
end
