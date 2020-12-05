p = traceplots(chn)

# posterior mean graphlap
#θ̄gl = vec(mean(θgl[bi_gl,:], dims=1))
θ̄dir = mat2vec(mean(θdir[bi_dir]))

hcat(θ̄gl, θ̄dir)

# hcat(θ̄gl, θ̄dir, θ̄mle[2:end])
# θ̄gl - θ̄mle[2:end]

# write probs to csv files
p_dir = write_binprobs(θ̄dir,truedatagen,binx,biny,"Dirichlet",θcopula)
p_gl = write_binprobs(θ̄gl,truedatagen,binx,biny,"graphLaplacian",θcopula)
error_dir_p = norm(p_dir[!,:pest]-p_dir[!,:ptrue])
error_gl_p = norm(p_gl[!,:pest]-p_gl[!,:ptrue])
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

# for traceplots, use iterates_dir and iterates_gl
#
# dirout = [Any[iterates_dir[iter][i,j], iter, "[$i,$j]"] for iter in 1:ITERdir, i in 1:m, j in 1:n][:]
# dirout_df = DataFrame(w=extractind(θiterates,1), iterate=extractind(θiterates,2), binID=extractind(θiterates,3))
# CSV.write("./out/diroutdf.csv",dirout_df)
#
# τindex = size(iterates_gl)[2]  # extract index where τ is stored in output of Turing
# glout = DataFrame(iterates_gl[:,[1,2,3,10,11,τindex],1])
# glout[!,:iterate] = 1:ITER_GL
# names!(glout,Symbol.(["iterate","H1","H2","H3","H10","H11","tau"]))
# CSV.write("./out/gloutdf.csv",glout) # el
