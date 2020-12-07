θ̄dir = mat2vec(mean(θdir[bi]))
θ̄gl = [mean(x[bi]) for x ∈ eachcol(θsave)]

#θ̄gl = vec(mean(θgl[bi_gl,:], dims=1))  # posterior mean graphlap using Turing
#p = traceplots(chn) # for Turing output


p2 = Plots.scatter(t[ind_yknown],y[ind_yknown],
color=:lightblue,label="y observed", legend=:outertopright)
Plots.scatter!(p2, t[ind_yunknown], y[ind_yunknown],color=:pink,
 label="y not observed", title = "Data",  xlabel="censoring time", ylabel="mark")

# traceplots for pcn
tr1 = Plots.plot(θsave[:,1], label="θ[1]")
tr2 = Plots.plot(θsave[:,10], label="θ[10]")
tr3 = Plots.plot(θsave[:,11], label="θ[11]")
tr4 = Plots.plot(τsave, label="τ")
lay = @layout [a b; c d]
Plots.plot(tr1, tr2, tr3, tr4, layout=lay)


# true binprobs
θ0, xx, yy = θtrue(truedatagen,binx,biny)
d = DataFrame(ptrue = θ0, Dirichlet =θ̄dir, graphLaplacian = θ̄gl,  x=xx, y=yy)
CSV.write("./out/binprobs.csv",d)
d

# compute Wasserstein distances
th0 = θ0; thdir = θ̄dir; thgl = θ̄gl
@rput th0 thdir thgl
R"""
wdir = wasserstein1d(th0, thdir)
wgl =  wasserstein1d(th0, thgl)
"""
@rget wdir wgl
println(wdir/wgl)

#----------------------------------------------------------------------------------------------
# write observations to csv file
yobserved = fill("yes",nsample)
yobserved[ind_yunknown] .= "no"
d = DataFrame(x=x,y=y,t=t,yobserved=yobserved)
CSV.write("./out/observations.csv",d)

#----------------------------------------------------------------------------------------------
# write info to file
facc = open("./out/info.txt","w")
    write(facc, "Data choice: ",truedatagen,"\n")
    write(facc, "Sample size: ", string(nsample), "\n")
    write(facc, "Fraction in data with mark observed: ",string(length(ind_yknown)/nsample),"\n\n")

    write(facc, "Hyperpar Dirichlet prior: ", string(ps), "\n\n")

    write(facc, "Sampler: ", string(sp),"\n")
    write(facc, "Number of iterations for Dirichlet prior: ",string(ITERdir),"\n")
    write(facc, "Number of burnin iterations for Dirichlet prior: ",string(BIdir),"\n")

    write(facc, "Number of iterations for graphLaplacian prior: ",string(ITERgl),"\n")
    write(facc, "Number of burnin iterations for graphLaplacian prior: ",string(BIgl),"\n")


    write(facc, "RootSquareError Dirichlet binprobs: ", string(error_dir_p),"\n")
    write(facc, "RootSquareError Graph Laplacian binprobs: ", string(error_gl_p),"\n\n")


close(facc)
