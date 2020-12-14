using JLD2

workdir = @__DIR__
cd(workdir)

nsample = 500#100



# example 1 (simple)
dist = Xplusy()

#dist = Uniform2D()



# example 2 (more interesting)
dist = Mixture(0.7,X2plusyRev(), X2plusy())


x, y, t, ind_yknown, ind_yunknown = gencensdata(dist, nsample)
JLD2.@save "data1.jld2" dist nsample x y t ind_yknown ind_yunknown
