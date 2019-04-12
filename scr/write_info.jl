minl = -0.25
maxl = 0.25
error_dir_p = plotting_p(θpostmean_dir,truedatagen,binx,biny,"Dirichlet";mincol_lim=minl,maxcol_lim=maxl)
error_gl_p = plotting_p(θpostmean_gl,truedatagen,binx,biny,"graphLaplacian";mincol_lim=minl,maxcol_lim=maxl)

minl = -2
maxl = 2
error_dir_d = plotting_d(θpostmean_dir/binarea,truedatagen,binx,biny, "Dirichlet";mincol_lim=minl,maxcol_lim=maxl)
error_gl_d = plotting_d(θpostmean_gl/binarea,truedatagen,binx,biny, "graphLaplacian";mincol_lim=minl,maxcol_lim=maxl)


yobserved = fill("yes",N)
yobserved[ind_yunknown] .= "no"
d = DataFrame(x=x,y=y,t=t,yobserved=yobserved)

facc = open("./out/info.txt","w")
    write(facc, "Data choice: ",truedatagen,"\n")
    write(facc, "Sample size: ", string(N), "\n")
    write(facc, "Fraction in data with mark observed: ",string(length(ind_yknown)/N),"\n\n")

    write(facc, "Hyperpar Dirichlet prior: ", string(priorscale), "\n\n")

    write(facc, "Sampler: ", string(gg),"\n")
    write(facc, "Number of iterations: ",string(ITER),"\n")
    write(facc, "Number of burnin iterations: ",string(BI),"\n")

    write(facc, "RootSquareError Dirichlet binprobs: ", string(error_dir_p),"\n")
    write(facc, "RootSquareError Graph Laplacian binprobs: ", string(error_gl_p),"\n\n")

    write(facc, "RootSquareError Dirichlet density: ", string(error_dir_d),"\n")
    write(facc, "RootSquareError Graph Laplacian density: ", string(error_gl_d),"\n")

close(facc)

# write observations to csv file
CSV.write("./out/observations.csv",d)
