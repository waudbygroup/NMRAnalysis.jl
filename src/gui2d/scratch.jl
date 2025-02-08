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

gui!(expt)

##
expt = HetNOEExperiment([
        "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/21/pdata/231",
        "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/21/pdata/232",
        "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/23/pdata/231",
        "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/23/pdata/232",
        "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/25/pdata/231",
        "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/25/pdata/232",
    ], [false, true, false, true, false, true])
gui!(expt)


##
expt = PREExperiment([
        "/Users/chris/NMR/crick-800/chris_hewl_231106/18/",
        "/Users/chris/NMR/crick-800/chris_hewl_231106/44/",
        "/Users/chris/NMR/crick-800/chris_hewl_231106/34/",
        "/Users/chris/NMR/crick-800/chris_hewl_231106/24/"
    ],
    [0., 1., 2., 5.],
    :hmqc,
    0.015)

gui!(expt)
