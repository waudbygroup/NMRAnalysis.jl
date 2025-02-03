using Pkg
pkg"activate ."
using NMRAnalysis
using GLMakie

ENV["JULIA_DEBUG"]=NMRAnalysis

##
expt = NMRAnalysis.GUI2D.RelaxationExperiment("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11",
    "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/lists/vd/t2-vd.cw")
state = NMRAnalysis.GUI2D.preparestate(expt)
# contour(expt.specdata.x[1], expt.specdata.y[1], expt.specdata.z[1])

NMRAnalysis.GUI2D.addpeak!(expt, Point2f(7.99,123.1), "X1")
NMRAnalysis.GUI2D.addpeak!(expt, Point2f(8.06,122.75), "X2")
# NMRAnalysis.GUI2D.addpeak!(expt, Point2f(7.98,122.45), "X3")
# NMRAnalysis.GUI2D.addpeak!(expt, Point2f(8.01,122.45), "X4")

expt.isfitting[] = true

NMRAnalysis.GUI2D.gui!(state, expt)

##
NMRAnalysis.GUI2D.deletepeak!(expt, 4)
