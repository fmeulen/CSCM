using JLD2

workdir = @__DIR__
cd(workdir)

Random.seed!(5)
nsample = 200

#dist = Uniform2D()

# example 1 (simple)
dist = Xplusy()
x, y, t, ind_yknown, ind_yunknown = gencensdata(dist, nsample)
JLD2.@save "data_Xplusy.jld2" dist nsample x y t ind_yknown ind_yunknown


# example 2 (more interesting)
dist = Mixture(0.7,X2plusyRev(), X2plusy())
x, y, t, ind_yknown, ind_yunknown = gencensdata(dist, nsample)
JLD2.@save "data_Mixture.jld2" dist nsample x y t ind_yknown ind_yunknown
