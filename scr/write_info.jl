#θ̄dir = mat2vec(mean(θdir[bi]))
θ̄dir = [mean(x[bi]) for x ∈ eachcol(θdir)]
θ̄gl = [mean(x[bi]) for x ∈ eachcol(θgl)]

# traceplots for pcn
# tr1 = Plots.plot(θsave[:,1], label="θ[1]")
# tr2 = Plots.plot(θsave[:,10], label="θ[10]")
# tr3 = Plots.plot(θsave[:,11], label="θ[11]")
# tr4 = Plots.plot(log.(τsave), label="log(τ)")
# lay = @layout [a b; c d]
# P = Plots.plot(tr1, tr2, tr3, tr4, layout=lay)
#savefig(P, "./out/traceplots_pcn.pdf")

dtrace = DataFrame(iterate=1:IT,theta1=θgl[:,1], theta10=θgl[:,10], theta11=θgl[:,11], logtau=log.(τgl))
CSV.write("./out/tracepcn.csv",dtrace)


# true binprobs
θ0, xx, yy = binprob(dist,bins)
d = DataFrame(ptrue = θ0, Dirichlet =θ̄dir, graphLaplacian = θ̄gl,  x=xx, y=yy)
CSV.write("./out/binprobs.csv",d)

distdir = norm(θ0 - θ̄dir,1)
distgl = norm(θ0 - θ̄gl,1)
@show distgl/distdir

# write observations to csv file
yobserved = fill("yes",nsample)
yobserved[ind_yunknown] .= "no"
d = DataFrame(x=x,y=y,t=t,yobserved=yobserved)
CSV.write("./out/observations.csv",d)

# heatmaps
m, n = bins.m, bins.n
heatmap(vec2mat(θ̄dir,m,n))
heatmap(vec2mat(θ̄gl,m,n))
heatmap(vec2mat(θ̄dir-θ0,m,n))
heatmap(vec2mat(θ̄gl-θ0,m,n))

# write info to file
facc = open("./out/info.txt","w")
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

    write(facc, "Fraction of accepted τ-update steps: ", string(accgl[2]/IT),"\n")
close(facc)
