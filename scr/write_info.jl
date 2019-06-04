######################## Trace plots for gl prior #######################
trace_df = DataFrame(th1=θitgl[:,1],th2=θitgl[:,2],th3=θitgl[:,3],
           th10=θitgl[:,10], th11=θitgl[:,11], tau = τgl )
nr, nc = size(trace_df)
@rput trace_df
@rput nr; @rput nc
R"""
library(tidyverse)
print(colnames(trace_df))
d <- trace_df %>% gather(value=value, key=coefficient) %>% mutate(iterate=rep(1:nr,nc))
d$coefficient <- factor(d$coefficient, labels=c("theta[1]","theta[2]","theta[3]",
                    "theta[10]","theta[11]","tau"))
d %>% ggplot(aes(x=iterate,y=value)) + geom_path() +
    facet_wrap(~coefficient,labeller=label_parsed,scales='free') +
    ylab("") + theme_light()+
    theme(strip.text.x = element_text(size=12,angle=0))
        ggsave("./out/traceplots.pdf",width=8, height=5)
"""

####################### binprob errors plots #######################
df_error_p = prep_plotting_p(θpm_dir,θpm_gl,truedatagen,binx,biny,"",θ)
error_p = plotting_p(df_error_p)

####################### density errors plots #######################
df_error_d = prep_plotting_d(θpm_dir/binarea,θpm_gl/binarea,truedatagen,binx,biny,θ ;gridN=200)
error_d = plotting_d(df_error_d)

# plot_truedensity(df_error_gl_d)

######################## write observations to csv file #######################
yobserved = fill("yes",N)
yobserved[ind_yunknown] .= "no"
dobs = DataFrame(x=x,y=y,t=t,yobserved=yobserved)
CSV.write("./out/observations.csv",dobs)

####################### write info to file #######################
facc = open("./out/info.txt","w")
    write(facc, "Data choice: ",truedatagen,"\n")
    write(facc, "Sample size: ", string(N), "\n")
    write(facc, "Fraction in data with mark observed: ",string(length(ind_yknown)/N),"\n\n")

    write(facc, "Hyperpar Dirichlet prior: ", string(priorscale), "\n\n")
    write(facc, "Number of iterations Dir: ",string(ITER),"\n")
    write(facc, "Number of burnin iterations Dir: ",string(BI),"\n")

    write(facc, "Sampler for graphLap: ", string(gg),"\n")
    write(facc, "Number of iterations graphLap: ",string(ITER_NUTS),"\n")
    write(facc, "Number of burnin iterations graphLap: ",string(BI_NUTS),"\n")

    write(facc, "RootSquareError Dirichlet binprobs: ", string(error_p[:error][1]),"\n")
    write(facc, "RootSquareError Graph Laplacian binprobs: ", string(error_p[:error][2]),"\n\n")

    write(facc, "RootSquareError Dirichlet density: ", string(error_d[:error][1]),"\n")
    write(facc, "RootSquareError Graph Laplacian density: ", string(error_d[:error][2]),"\n")

close(facc)
