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
close(facc)

# write observations to csv file
CSV.write("./out/observations.csv",d)
