# write probs to csv files
p_dir = write_binprobs(θpostmean_dir,truedatagen,binx,biny,"Dirichlet",θcopula)
p_gl = write_binprobs(θpostmean_gl,truedatagen,binx,biny,"graphLaplacian",θcopula)
error_dir_p = norm(p_dir[:pest]-p_dir[:ptrue])
error_gl_p = norm(p_gl[:pest]-p_gl[:ptrue])
#----------------------------------------------------------------------------------------------
# write observations to csv file
yobserved = fill("yes",NSAMPLE)
yobserved[ind_yunknown] .= "no"
d = DataFrame(x=x,y=y,t=t,yobserved=yobserved)
CSV.write("./out/observations.csv",d)

#----------------------------------------------------------------------------------------------
# write info to file
facc = open("./out/info.txt","w")
    write(facc, "Data choice: ",truedatagen,"\n")
    write(facc, "Sample size: ", string(NSAMPLE), "\n")
    write(facc, "Fraction in data with mark observed: ",string(length(ind_yknown)/NSAMPLE),"\n\n")

    write(facc, "Hyperpar Dirichlet prior: ", string(ps), "\n\n")

    write(facc, "Sampler: ", string(sp),"\n")
    write(facc, "Number of iterations for Dirichlet prior: ",string(ITER_DIR),"\n")
    write(facc, "Number of burnin iterations for Dirichlet prior: ",string(BI_DIR),"\n")

    write(facc, "Number of iterations for graphLaplacian prior: ",string(ITER_GL),"\n")
    write(facc, "Number of burnin iterations for graphLaplacian prior: ",string(BI_GL),"\n")


    write(facc, "RootSquareError Dirichlet binprobs: ", string(error_dir_p),"\n")
    write(facc, "RootSquareError Graph Laplacian binprobs: ", string(error_gl_p),"\n\n")


close(facc)
