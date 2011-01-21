import testcr

def mtr():
    nruns = 3
    ntrajlist = [1,2]
    nchannels = 256
    npoints = 21
    tcr = testcr.mk_tcr()
    seed = 0
    saveList = []
    for nt in ntrajlist:
        saveList.append([])
        for nr in range(nruns):
            seed +=1
            tcr.datagen.fill(nchannels,seed,nt)
            tcr.generator.fitnesslist.o(0).npoints(npoints)
            tcr.mrf.opt.set_optimizer("BFGSWrap")
            tcr.compute(nchannels, seed, nt, run=4)
            saveList[-1].append(tcr.save())
    return saveList
