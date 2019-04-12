# prepare df
df_error_dir_p = prep_plotting_p(θpostmean_dir,truedatagen,binx,biny,"Dirichlet")
df_error_gl_p = prep_plotting_p(θpostmean_gl,truedatagen,binx,biny,"graphLaplacian")
# set colour limits
minl, maxl = extrema(  hcat(df_error_dir_p.pest-df_error_dir_p.ptrue,df_error_gl_p.pest-df_error_gl_p.ptrue))
# make plots
error_dir_p = plotting_p(df_error_dir_p,"Dirichlet";mincol_lim=minl,maxcol_lim=maxl)
error_gl_p = plotting_p(df_error_gl_p,"graphLaplacian";mincol_lim=minl,maxcol_lim=maxl)

# similarly for densities
df_error_dir_d = prep_plotting_d(θpostmean_dir/binarea,truedatagen,binx,biny, "Dirichlet")
df_error_gl_d = prep_plotting_d(θpostmean_gl/binarea,truedatagen,binx,biny, "graphLaplacian")
minl, maxl = extrema(  hcat(df_error_dir_d.dest-df_error_dir_d.dtrue,df_error_gl_d.dest-df_error_gl_d.dtrue))
error_dir_d = plotting_d(df_error_dir_d, "Dirichlet";mincol_lim=minl,maxcol_lim=maxl)
error_gl_d = plotting_d(df_error_gl_d, "graphLaplacian";mincol_lim=minl,maxcol_lim=maxl)

# write observations to csv file
yobserved = fill("yes",N)
yobserved[ind_yunknown] .= "no"
d = DataFrame(x=x,y=y,t=t,yobserved=yobserved)
CSV.write("./out/observations.csv",d)


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
