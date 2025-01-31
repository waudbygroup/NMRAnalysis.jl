using Pkg
pkg"activate ."
using NMRAnalysis
using GLMakie

ENV["JULIA_DEBUG"]=NMRAnalysis

##
expt = NMRAnalysis.GUI2D.RelaxationExperiment("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11",
    "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/lists/vd/t2-vd.cw")

# contour(expt.specdata.x[1], expt.specdata.y[1], expt.specdata.z[1])

expt.isfitting[] = true
##

NMRAnalysis.GUI2D.addpeak!(expt, Point2f(7.99,123.1), "X1")
NMRAnalysis.GUI2D.addpeak!(expt, Point2f(8.06,122.75), "X2")
NMRAnalysis.GUI2D.addpeak!(expt, Point2f(7.98,122.45), "X3")
NMRAnalysis.GUI2D.addpeak!(expt, Point2f(8.01,122.45), "X4")
##
NMRAnalysis.GUI2D.deletepeak!(expt, 4)
##
f = Figure()
ax1 = Axis(f[1,1])
ax2 = Axis(f[1,2])
ax3 = Axis(f[1,3])
linkaxes!(ax1, ax2, ax3)
contour!(ax1, expt.specdata.x[1], expt.specdata.y[1], expt.specdata.z[1])
contour!(ax2, expt.specdata.x[1], expt.specdata.y[1], expt.specdata.z[1], color=:blue, levels=20 .*[5,10,15,20,25,30])
contour!(ax2, expt.specdata.x[1], expt.specdata.y[1], expt.specdata.zfit[][1], color=:red, levels=20 .*[5,10,15,20,25,30])
heatmap!(ax3,expt.specdata.x[1], expt.specdata.y[1], expt.specdata.mask[][1])

f