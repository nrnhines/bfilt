import neuron
import neuron.gui
import init
from neuron import h
import noisegen
import scipy.stats as stats
import numpy

class TestCR(object):
    def __init__(self,n,seed,modelses,datagenhoc):
        self.alpha = 0.05
        self.n = n
        h.load_file(modelses)
        self.N = h.List("PythonObject").o(0)
        self.true = self.N.getParm()
        self.modelses = modelses
        if n == 0:
          G = noisegen.Gen(self.N)
          G.reseed(seed)
          self.seed = seed
          self.Data = G.datasim()
        else:
          h.load_file(datagenhoc)
          tvec = h.Vector(self.N.Eve.collectionTimes)
          vec = h.ch4ssdata(n, seed, tvec, self.true)
          self.Data = []
          for i in range(len(vec)):
            self.Data.append(numpy.matrix(vec[i]))
        h.topology()
        ss = h.Vector()
        h.cvode.states(ss)
        ss.printf()
        self.N.overwrite(self.Data)
        self.tl = self.N.likelihood()
        print self.tl
        h.attr_praxis(seed)
        h.MulRunFitter[0].efun()
        self.mle = self.N.getParm()
        self.ml = self.N.likelihood()
        self.CS = 2.0*(self.tl - self.ml)
        self.pValue = stats.chisqprob(self.CS,self.true.size())
        self.covers = (self.pValue >= self.alpha)
        self.N.setParm(self.true)

def run(nruns=1,nchannels=100,modelses="ch4.ses",datagenhoc="ch4ssdatagen.hoc"):
    TCRs = []
    for i in range(nruns):
        TCRi = TestCR(nchannels,i+1,modelses,datagenhoc)
        TCRs.append(TCRi)
        if TCRi.covers:
            print i, 'IN:', TCRi.pValue
        else:
            print i, 'OUT:', TCRi.pValue
    return TCRs
