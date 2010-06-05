import neuron
import neuron.gui
import init
from neuron import h
import noisegen
import scipy.stats as stats

class TestCR(object):
    def __init__(n,ses,seed):
        assert(n == 0)
        self.alpha = 0.05
        self.n = n
        h.load_file(ses)
        self.ses = ses
        self.N = h.List("PythonObject").o(0)
        self.G = noisegen.Gen(self.N)
        self.G.reseed(seed)
        self.seed = seed
        self.Data = self.G.datasim()
        self.N.overwrite(self.Data)
        self.true = self.N.getParm()
        self.tl = self.N.likelihood()
        h.MulRunFitter[0].efun()
        self.mle = self.N.getParm()
        self.ml = self.N.likelihood()
        self.CS = 2.0*(self.tl - self.ml)
        self.pValue = stats.chisqprob(self.CS,self.true.size())
        self.covers = (self.pValue >= self.alpha)