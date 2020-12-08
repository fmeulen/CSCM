θ̄dir = mat2vec(mean(θdir[bi]))
θ̄gl = [mean(x[bi]) for x ∈ eachcol(θsave)]

#θ̄gl = vec(mean(θgl[bi_gl,:], dims=1))  # posterior mean graphlap using Turing
#p = traceplots(chn) # for Turing output


p2 = Plots.scatter(t[ind_yknown],y[ind_yknown],
color=:lightblue,label="y observed", legend=:outertopright)
Pobs = Plots.scatter!(p2, t[ind_yunknown], y[ind_yunknown],color=:pink,
 label="y not observed", title = "Data",  xlabel="censoring time", ylabel="mark")
savefig(Pobs, "./out/observations.pdf")

# traceplots for pcn
tr1 = Plots.plot(θsave[:,1], label="θ[1]")
tr2 = Plots.plot(θsave[:,10], label="θ[10]")
tr3 = Plots.plot(θsave[:,11], label="θ[11]")
tr4 = Plots.plot(log.(τsave), label="log(τ)")
lay = @layout [a b; c d]
P = Plots.plot(tr1, tr2, tr3, tr4, layout=lay)

savefig(P, "./out/traceplots_pcn.pdf")



# true binprobs
θ0, xx, yy = θtrue(truedatagen,binx,biny)
d = DataFrame(ptrue = θ0, Dirichlet =θ̄dir, graphLaplacian = θ̄gl,  x=xx, y=yy)
CSV.write("./out/binprobs.csv",d)

# compute Wasserstein distances
th0 = θ0; thdir = θ̄dir; thgl = θ̄gl
@rput th0 thdir thgl
R"""
library(transport)
wdir = wasserstein1d(th0, thdir)
wgl =  wasserstein1d(th0, thgl)
"""
@rget wdir wgl
@show ratio = round(wgl/wdir; digits=3)

# write observations to csv file
yobserved = fill("yes",nsample)
yobserved[ind_yunknown] .= "no"
d = DataFrame(x=x,y=y,t=t,yobserved=yobserved)
CSV.write("./out/observations.csv",d)

# heatmaps
heatmap(mean(θdir[bi]))
heatmap(vec2mat(θ̄gl,m,n))
heatmap(mean(θdir[bi]) - vec2mat(θ0,m,n))
heatmap(vec2mat(θ̄gl-θ0,m,n))

# write info to file
facc = open("./out/info.txt","w")
    write(facc, "Data choice: ",truedatagen,"\n")
    write(facc, "Sample size: ", string(nsample), "\n")
    write(facc, "Fraction in data with mark observed: ",string(length(ind_yknown)/nsample),"\n\n")

    write(facc, "Hyperpar Dirichlet prior: ", string(ps), "\n\n")

    write(facc, "Number of iterations: ",string(IT),"\n")
    write(facc, "Number of burnin iterations: ",string(BI),"\n")

    if false
        write(facc, "Turing Sampler: ", string(sp),"\n")
        write(facc, "Number of iterations for Turing: ",string(ITERgl),"\n")
        write(facc, "Number of burnin iterations for Turing: ",string(BIgl),"\n")
    end

    write(facc, "Wasserstein distance binprobs Dirichlet: ", string(wdir),"\n")
    write(facc, "Wasserstein distance binprobs Graph Laplacian: ", string(wgl),"\n")
    write(facc, "Ratio: ", string(ratio), "\n\n")

    write(facc, "Parameter rho for pCN:", string(ρ), "\n")
    write(facc, "Fraction of accepted pCN steps: ", string(acc[1]/IT),"\n\n")

    write(facc, "Fraction of accepted τ-update steps: ", string(acc[2]/IT),"\n")
close(facc)
