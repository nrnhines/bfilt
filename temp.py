from testcr import *
tcr = mk_tcr()
tcr.datagen.fill(50, 1, 3)
tcr.generator.fitnesslist.o(0).npoints(21)
cvodewrap.fs.use_fixed_step = 1.0
