import testcr

def mtr():
    nruns = 3
    ntrajlist = [1,2]
    nchannels = 256
    npoints = 21
    tcr = testcr.mk_tcr()
    seed = 0
    pValues = []
    for nt in ntrajlist:
        pValues.append([])
        for nr in range(nruns):
            seed +=1
            tcr.datagen.fill(nchannels,seed,nt)
            tcr.generator.fitnesslist.o(0).npoints(npoints)
            tcr.compute(self, nchannels, seed, nt, run=4)
            pValues[-1].append(tcr.APFpValue)
