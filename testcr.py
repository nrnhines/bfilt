import neuron
#import neuron.gui
import init
from neuron import h
import noisegen
import scipy.stats as stats
import numpy
import numdifftools as nd
import math
import pickle
import svd
import cvodewrap
import fitglobals

def first(modelses):
    h.load_file(modelses)
    Z = h.MulRunFitter[0].p.pf.parmlist
    #Z.append(h.RunFitParm("nb.Eve.Sto.scale",1,1e-9,1e9,1,1))

class TestCR(object):
    def __init__(self,n,seed,modelses,datagenhoc,run=2):
        self.alpha = 0.05
        self.n = n  # number of channels
        h.load_file(modelses)
        mrflist = h.List("MulRunFitter")
        mrf = mrflist.o(int(mrflist.count())-1)
        self.N = mrf.p.pf.generatorlist.o(0).gen.po
        h('objref nb')
        h.nb = mrf.p.pf.generatorlist.o(0).gen.po
        cvodewrap.fs.panel()
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
          vec = h.ch3ssdata(n, seed, tvec, self.true, self.N.rf.fitnesslist.o(0))
          self.Data = []
          for i in range(len(vec)):
            self.Data.append(numpy.matrix(vec[i]))
        if fitglobals.verbose: h.topology()
        ss = h.Vector()
        cvodewrap.states(ss)
        if fitglobals.verbose: ss.printf()
        self.N.overwrite(self.Data)
        # self.tl = self.N.likelihood()
        # print self.tl
        self.Z = h.MulRunFitter[0].p.pf.parmlist
        if fitglobals.verbose: print "ASSUMES PARAMETERS 0,1 main parameters rest NUISANCE"
        # DESTROY
        while self.Z.count() > 2:
            self.Z.remove(self.Z.count()-1)
        h.FitnessGenerator1[0].use_likelihood=0
        self.pre0Parm = self.N.getParmVal()
        h.MulRunFitter[0].efun()
        self.post0Parm = self.N.getParmVal()
        h.FitnessGenerator1[0].use_likelihood=1
        foo = h.RunFitParm("nb.Eve.Sto.scale")
        foo.set("nb.Eve.Sto.scale",1,1e-9,1e9,1,1)
        self.Z.append(foo)
        self.Z.o(0).doarg = 0
        self.Z.o(1).doarg = 0
        self.Z.o(2).doarg = 1
        h.attr_praxis(seed)
        #print 'SIZE =', self.N.getParm().size()
        if run == 0:
          return
        self.pre1Parm = self.N.getParmVal()
        h.MulRunFitter[0].efun()
        self.post1Parm = self.N.getParmVal()
        self.otle = self.N.getParmVal()
        self.otml = self.N.likelihood()  #optimized true maximum likelihood
        if run == 1:
          return
        self.Z.o(0).doarg = 1
        self.Z.o(1).doarg = 1
        self.Z.o(2).doarg = 1
        h.attr_praxis(seed)
        self.pre2Parm = self.N.getParmVal()
        h.MulRunFitter[0].efun()
        self.post2Parm = self.N.getParmVal()
        self.mle = self.N.getParmVal()
        self.ml = self.N.likelihood()
        self.CS = 2.0*(self.otml - self.ml)
        self.pValue = stats.chisqprob(self.CS,self.true.size())
        self.covers = (self.pValue >= self.alpha)
        if run == 2:
          return
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

class TestRao(TestCR):
    def __init__(self,n,seed,modelses,datagenhoc):
        self.n = n
        h.load_file(modelses)
        mrflist = h.List("MulRunFitter")
        mrf = mrflist.o(int(mrflist.count())-1)
        self.N = mrf.p.pf.generatorlist.o(0).gen.po
        h('objref nb')
        h.nb = mrf.p.pf.generatorlist.o(0).gen.po
        cvodewrap.fs.panel()
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
            vec = h.ch3ssdata(n, seed, tvec, self.true, self.N.rf.fitnesslist.o(0))
            self.Data = []
            for i in range(len(vec)):
                self.Data.append(numpy.matrix(vec[i]))
        if fitglobals.verbose: h.topology()
        ss = h.Vector()
        cvodewrap.states(ss)
        if fitglobals.verbose: ss.printf()
        self.N.overwrite(self.Data)
        self.H = numpy.matrix(self.Hessian())

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

def start(seed=1, nchannels=50, modelses="ch3_101p.ses", datagenhoc="ch3ssdatagen.hoc"):
    return TestCR(nchannels,seed,modelses,datagenhoc,run=0)

def prin():
    v = h.Vector()
    x = h.pval_praxis(0, v)
    n = len(v)
    pval = [x]
    paxis = [numpy.array(v)]
    for i in range(1, n):
        pval.append(h.pval_praxis(i, v))
        paxis.append(numpy.array(v))
    return (pval, paxis)

def onerun(seed=1, nchannels=50, modelses="ch3_101p.ses", datagenhoc="ch3ssdatagen.hoc"):
    cvodewrap.fs.use_fixed_step = 1.0
    tt = h.startsw()
    pc = h.ParallelContext()
    id = int(pc.id())
    print "%d start seed=%d nchannels=%d"%(id, seed, nchannels)
    r = TestCR(nchannels,seed,modelses,datagenhoc, run=2)
    tt = h.startsw() - tt
    print "%d finish walltime=%g seed=%d nchannels=%d"%(id, tt, seed, nchannels)
    #return value suitable for bulletin board
    return (tt, (seed, nchannels), numpy.array(r.otle), r.otml, numpy.array(r.mle), r.ml, prin())


def runRao(nruns=1,nchannels=50,modelses="ch3_101p.ses",datagenhoc="ch3ssdatagen.hoc"):
    TCRs = []
    for i in range(nruns):
        TCRi = TestRao(nchannels,i+1,modelses,datagenhoc)
        TCRs.append(TCRi)
        print TCRi.H
    return TCRs

def calcRaoCov(TCRs):
    EH = TCRs[0].H*0.0
    w = 1.0/len(TCRs)
    for T in TCRs:
        EH += w*T.H
    return EH.I

def run(nruns=1,nchannels=50,modelses="ch3_101p.ses",datagenhoc="ch3ssdatagen.hoc"):
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

def parrun(nruns=0,nchannels=50,modelses="ch3_101p.ses",datagenhoc="ch3ssdatagen.hoc"):
    pc = h.ParallelContext()
    pc.runworker()
    if nruns == 0:
        nruns = int(pc.nhost())
    for i in range(nruns):
        pc.submit(onerun, i+1, nchannels, modelses, datagenhoc)
    f = open('results_101p_'+str(nchannels)+'.'+str(nruns), "w")
    pickle.dump(nruns, f); f.flush()
    i = 0
    while (pc.working() != 0.0):
        r = pc.pyret()
        pickle.dump(r, f); f.flush()
        i += 1
        print '%d seed %d'%(i,r[1][0])
    f.close()
    pc.done()
    h.quit()

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
    first('ch3_101p.ses')
    T = run(3,10000,'ch3_101p.ses')
    pickelWOH(T,'TCRtry.pkl')
