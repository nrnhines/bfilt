import testcr
import pickle

def readmtr(filename):
    f = open(filename,'r')
    mtr = pickle.load(f)
    f.close()
    return mtr

class MTR(object):
    def __init__(self, seed_offset=0):
        self.seed_offset = seed_offset
        self.nruns = 100
        self.ntrajlist = [1,4,16]
        self.nchannels = 50
        self.npoints = 21
        self.output = None

    def run(self):
        if not self.output == None:
            return
        tcr = testcr.mk_tcr()
        seed = 0
        self.output = []
        for nt in self.ntrajlist:
            self.output.append([])
            for nr in range(self.nruns):
                seed +=1
                tcr.datagen.fill(self.nchannels,seed+self.seed_offset,nt)
                tcr.generator.fitnesslist.o(0).npoints(self.npoints)
                tcr.mrf.opt.set_optimizer("BFGSWrap")
                tcr.compute(self.nchannels,seed+self.seed_offset, nt, run=4)
                self.output[-1].append(tcr.save())

    def save(self,filename):
        f = open(filename,'w')
        pickle.dump(self,f)
        f.close()

