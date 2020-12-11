using JLD2

workdir = @__DIR__
cd(workdir)

nsample = 100


dist = X2plusy()
x, y, t, ind_yknown, ind_yunknown = gencensdata(dist, nsample)


JLD2.@save "data1.jld2" dist nsample x y t ind_yknown ind_yunknown





#dist = Mixture(0.3, XplusyRev(), Xplusy())
#dist = Uniform2D()
