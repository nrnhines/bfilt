import neuron
import neuron.gui
import init
from neuron import h
import noisegen
import scipy.stats as stats
import numpy
import numdifftools as nd
import math
import pickle
import svd

class TestCR(object):
    def __init__(self,n,seed,modelses,datagenhoc):
        self.alpha = 0.05
        self.n = n
        h.load_file(modelses)
        self.N = h.List("PythonObject").o(0)
        self.true = self.N.getParm()
        self.modelses = modelses
        self.seed = seed
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
        self.H = numpy.matrix(self.Hessian())
        svdList = svd.svd(numpy.array(self.H))[1]
        self.precision = 0.0
        for sl in svdList:
            sl_positive = max(sl,1e-14)
            self.precision += math.log(sl_positive)
        self.N.setParm(self.true)
        self.likefailed = self.N.likefailed
        self.N.likefailed = False

    def evalFun(self, p):
        saveParm = self.N.getParm()
        newParm = self.N.getParm()
        for i in range(len(p)):
            newParm.x[i] = p[i]
        self.N.setParm(newParm)
        L = self.N.likelihood()
        self.N.setParm(saveParm)
        return L

    def Hessian(self):
        # Read current values of parameters (hopefully MLEs)
        parm = self.N.getParm()
        parmList = []
        for i in range(int(parm.size())):
            parmList.append(parm[i])
        # Create (likelihood) inline function
        LamFun = lambda p: self.evalFun(p)
        # Create Hessian (of likelihood) inline function
        HessFun = nd.Hessian(LamFun)
        # Evaluate Hessian and return
        return numpy.matrix(HessFun(parmList))

class WOHoc(object):
    def __init__(self, WH):
        self.alpha = WH.alpha
        self.n = WH.n
        self.true = numpy.matrix(WH.true)
        self.modelses = WH.modelses
        self.seed = WH.seed
        self.tl =WH.tl
        self.mle = numpy.matrix(WH.mle)
        self.ml = WH.ml
        self.CS = WH.CS
        self.pValue = WH.pValue
        self.covers = WH.covers
        self.H = WH.H
        self.precision =WH.precision
        self.likefailed = WH.likefailed

def onerun(seed, nchannels, modelses, datagen):
    return TestCR(nchannels,seed,modelses,datagenhoc)

def run(nruns=1,nchannels=100,modelses="ch4.ses",datagenhoc="ch4ssdatagen.hoc"):
    global pc
    TCRs = []
    ncovers = 0
    for i in range(nruns):
        pc.submit(onerun, nchannels, i+1, modelses, datagenhoc)
    while (pc.working()):
        TCRi = pc.pyret()
        TCRs.append(TCRi)
        if TCRi.covers:
            print i, 'IN:', TCRi.pValue
            ncovers += 1
        else:
            print i, 'OUT:', TCRi.pValue
    print ncovers, 'covers out of', nruns
    return TCRs

def pickelWOH(TCRs,filename):
    T = []
    for i in range(len(TCRs)):
        T.append(WOHoc(TCRs[i]))
    f = open(filename,"w")
    pickle.dump(T,f)
    f.close()

def load(filename):
    f = open(filename,'r')
    T = pickle.load(f)
    f.close()
    return T

def batch(nrun=1000):
    T = run(nrun,0,'rc64.ses')
    pickelWOH(T,'T64RCK.pkl')

pc = ParallelContext()
pc.runworker()
batch(4)
pc.done()
h.quit()

