using Pkg
pkg"activate ."
using NMRAnalysis

expt = NMRAnalysis.GUI2D.RelaxationExperiment("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11",
    "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/lists/vd/t2-vd.cw")
contour(expt.specdata.x[1], expt.specdata.y[1], expt.specdata.z[1])
NMRAnalysis.GUI2D.addpeak!(expt, Point2f(7.99,123.1), "X1")
