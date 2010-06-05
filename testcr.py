import neuron
import neuron.gui
import init
from neuron import h
import noisegen
import scipy.stats as stats

class TestCR(object):
    def __init__(self,n,ses,seed):
        assert(n == 0)
        self.alpha = 0.05
        self.n = n
        h.load_file(ses)
        self.ses = ses
        self.N = h.List("PythonObject").o(0)
        G = noisegen.Gen(self.N)
        G.reseed(seed)
        self.seed = seed
        self.Data = G.datasim()
        self.N.overwrite(self.Data)
        self.true = self.N.getParm()
        self.tl = self.N.likelihood()
        h.MulRunFitter[0].efun()
        self.mle = self.N.getParm()
        self.ml = self.N.likelihood()
        self.CS = 2.0*(self.tl - self.ml)
        self.pValue = stats.chisqprob(self.CS,self.true.size())
        self.covers = (self.pValue >= self.alpha)
        self.N.setParm(self.true)

def run(n=100):
    TCRs = []
    for i in range(n):
        TCRi = TestCR(0,"ch3.ses",i)
        TCRs.append(TCRi)
        if TCRi.covers:
            print i, 'IN:', TCRi.pValue
        else:
            print i, 'OUT:', TCRi.pValue
    return TCRs