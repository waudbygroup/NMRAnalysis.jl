using Pkg
pkg"activate ."
using NMRAnalysis
using GLMakie

ENV["JULIA_DEBUG"]=NMRAnalysis

##
expt = RelaxationExperiment("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11",
    "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/lists/vd/t2-vd.cw")
# expt = RelaxationExperiment("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/17",
#     "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/17/lists/vd/t1.cw")

NMRAnalysis.GUI2D.gui!(expt)

