# visualisaton of the posterior
error_dir_p,  error_gl_p, error_dir_d,  error_gl_d = plotting(θpostmean_dir,truedatagen,binx,biny,binarea,θcopula, NSAMPLE)

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

    write(facc, "RootSquareError Dirichlet density: ", string(error_dir_d),"\n")
    write(facc, "RootSquareError Graph Laplacian density: ", string(error_gl_d),"\n")

close(facc)

#----------------------------------------------------------------------------------------------
if make_traceplots
    ##### trace plots for graphlap
    trace_df = DataFrame(th1=θiterates_gl[:,1],th2=θiterates_gl[:,2],th10=θiterates_gl[:,10],
    th11=θiterates_gl[:,11])
    nr, nc = size(trace_df)
    @rput trace_df
    @rput nr; @rput nc
    R"""
    library(ggplot2);    library(tidyverse);      theme_set(theme_light())
    trace_df %>% gather(value=value, key=coefficient) %>% mutate(iterate=rep(1:nr,nc)) %>%
        ggplot(aes(x=iterate,y=value)) + geom_path() + facet_wrap(~coefficient,scales='free')
    """

    #### trace plots for dirichlet
    θiterates = [Any[iterates_dir[iter][i,j], iter, "[$i,$j]"] for iter in 1:ITER, i in 1:m, j in 1:n][:]
    iterates_df = DataFrame(w=extractind(θiterates,1), iterate=extractind(θiterates,2), binID=extractind(θiterates,3))
    @rput iterates_df
    R"""
    library(ggplot2);    library(tidyverse);      theme_set(theme_light())
    iterates_df %>%  ggplot() + geom_line(aes(x=iterate,y=w)) + facet_wrap(~binID,scales="free")
    """
end
