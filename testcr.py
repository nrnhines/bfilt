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

def first(modelses):
    h.load_file(modelses)
    Z = h.MulRunFitter[0].p.pf.parmlist
    #Z.append(h.RunFitParm("nb.Eve.Sto.scale",1,1e-9,1e9,1,1))

class TestCR(object):
    def __init__(self,n,seed,modelses,datagenhoc):
        self.alpha = 0.05
        self.n = n  # number of channels
        h.load_file(modelses)
        self.N = h.List("PythonObject").o(0)
        h('objref nb')
        h.nb = h.List("PythonObject").o(0)
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
          vec = h.ch3ssdata(n, seed, tvec, self.true)
          self.Data = []
          for i in range(len(vec)):
            self.Data.append(numpy.matrix(vec[i]))
        h.topology()
        ss = h.Vector()
        h.cvode.states(ss)
        ss.printf()
        self.N.overwrite(self.Data)
        # self.tl = self.N.likelihood()
        # print self.tl
        self.Z = h.MulRunFitter[0].p.pf.parmlist
        print "ASSUMES PARAMETERS 0,1 main parameters rest NUISANCE"
        self.Z.append(h.RunFitParm("nb.Eve.Sto.scale",1,1e-9,1e9,1,1))
        self.Z.o(0).doarg = 0
        self.Z.o(1).doarg = 0
        h.attr_praxis(seed)
        print 'SIZE =', self.N.getParm().size()
        h.MulRunFitter[0].efun()
        # return
        self.otle = self.N.getParm()
        self.otml = self.N.likelihood()  #optimized true maximum likelihood
        self.Z.o(0).doarg = 1
        self.Z.o(1).doarg = 1
        h.attr_praxis(seed)
        h.MulRunFitter[0].efun()
        self.mle = self.N.getParm()
        self.ml = self.N.likelihood()
        self.CS = 2.0*(self.otml - self.ml)
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
        self.otml = WH.otml
        self.mle = numpy.matrix(WH.mle)
        self.ml = WH.ml
        self.CS = WH.CS
        self.pValue = WH.pValue
        self.covers = WH.covers
        self.H = WH.H
        self.precision =WH.precision
        self.likefailed = WH.likefailed

def run(nruns=1,nchannels=50,modelses="ch3.ses",datagenhoc="ch3ssdatagen.hoc"):
    TCRs = []
    ncovers = 0
    for i in range(nruns):
        TCRi = TestCR(nchannels,i+1,modelses,datagenhoc)
        TCRs.append(TCRi)
        # return TCRs
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

def batch():
    first('ch3.ses')
    T = run(3,10000,'ch3.ses')
    pickelWOH(T,'TCRtry.pkl')
