import neuron
# import neuron.gui
import init
from neuron import h
h.ivoc_style("*xvalue_format", "%.15g")
import noisegen
import scipy.stats as stats
import numpy
import numdifftools as nd
import math
import pickle
import svd
import cvodewrap
import fitglobals
import repro

MPF = True

def first(modelses):
    h.load_file(modelses)
    Z = h.MulRunFitter[0].p.pf.parmlist
    #Z.append(h.RunFitParm("nb.Eve.Sto.scale",1,1e-9,1e9,1,1))

def MulRunFitHandle():
    mrflist = h.List("MulRunFitter")
    mrf = mrflist.o(int(mrflist.count())-1)
    return mrf

def NrnBFiltHandle(mrf):
    N = mrf.p.pf.generatorlist.o(0).gen.po
    h('objref nb')
    h.nb = N
    return N

def printSomeInfo():
    if fitglobals.verbose: h.topology()
    ss = h.Vector()
    cvodewrap.states(ss)
    if fitglobals.verbose: ss.printf()

class TestCR(object):
    def destroyNuisanceParms(self):  #Destroy all but nMain=2 Parms
        while self.parmlist.count() > self.nMain:
            self.parmlist.remove(self.parmlist.count()-1)
        self.nuisanceParms = []
        self.nNuisance = 0

    def addNuisanceParm(self, parmName):
        foo = h.RunFitParm(parmName)
        foo.set(parmName,1,1e-9,1e9,1,1)
        self.parmlist.append(foo)

    def usingArgs(self, fit_the_real_parms=True, fit_the_nuisance_parms=False):
        for i in range(self.nMain):
            self.parmlist.o(i).doarg = fit_the_real_parms
        for i in range(self.nMain, len(self.parmlist)):
            self.parmlist.o(i).doarg = fit_the_nuisance_parms
        self.mrf.p.pf.def_parmlist_use()

    def __init__(self, multiple_run_fitter, nrnbfilt, datagen):
        # datagen is a Repro from repro.py and fills the multiple run fitter fitness
        # function with reproducible data via datagen.fill(n, seed) in such a
        # way that the likelihood calculation uses the latest data

        # If available, code tests if the confidence region covers the true parms
        self.alpha = 0.05  # p < alpha outside confidence region, 0.5: 95% CRs
        self.mrf = multiple_run_fitter
        self.N = nrnbfilt
        self.datagen = datagen
        self.trueParm = h.true_parameters #from exper_channel.hoc set by exper_data.hoc
        self.saveParm = self.N.getParm()
        assert(len(self.trueParm) == len(self.saveParm))
        self.nMain = len(self.saveParm)
        if fitglobals.verbose: print "ASSUMES PARAMETERS 0...nMain-1 main parameters rest NUISANCE"
        self.nNuisance = 1
        self.parmlist = self.mrf.p.pf.parmlist
        self.ef = self.mrf.p.run   # ef is a function
        self.nuisanceParms = []
        self.generator = self.mrf.p.pf.generatorlist.o(0).gen
        self.usingArgs(True, False)
        self.N.setParmVal(self.trueParm)
        self.addNuisanceParm("nb.Eve.Sto.scale")

    def get_pValue(self, old, new, size):
        CS = 2.0*(old - new)
        pValue = stats.chisqprob(CS, size)
        return pValue

    def compute(self, nchan, seed, n_trajec, run=4):
        # The "run" parameter controls when the fitting stops
        # If run = 0 doesn't fit
        self.lastrun = run
        self.lastnchan = nchan
        self.lastseed = seed
        self.lastn_trajec = n_trajec
        self.datagen.fill(nchan, seed, n_trajec)
        cvodewrap.fs.panel()
        printSomeInfo()
        # self.tl = self.N.likelihood()
        # print self.tl
        if run == 0:
            return

        # Nuisance Fit At True
        self.usingArgs(False, True)
        h.attr_praxis(seed)
        #print 'SIZE =', self.N.getParm().size()
        self.preTrueParm = self.N.getParmVal()
        self.preTruef = self.ef()
        self.mrf.efun()
        self.postTrueParm = self.N.getParmVal()
        self.postTruef = self.ef()
        self.otle = self.N.getParmVal()
        self.otml = self.N.likelihood()  #optimized true maximum likelihood
        if run == 1:
          return

        # Square Norm Fit
        self.usingArgs(True, False)
        self.generator.use_likelihood=0
        h.attr_praxis(seed)
        self.preSNFParm = self.N.getParmVal()
        self.preSNFf = self.ef()
        self.mrf.efun()
        self.postSNFParm = self.N.getParmVal()
        self.postSNFf = self.ef()
        self.generator.use_likelihood=1
        if run == 2:
            return

        # Main Parm Fit:
        global MPF
        if MPF:
            h.attr_praxis(seed)
            self.preMPFParm = self.N.getParmVal()
            self.preMPFf = self.ef()
            self.mrf.efun()
            self.postMPFParm = self.N.getParmVal()
            self.postMPFf = self.ef()
            self.MPFpValue = self.get_pValue(self.postTruef, self.postMPFf, self.trueParm.size())
        if run == 3:
            return

        # All Parm Fit
        self.usingArgs(True, True)
        h.attr_praxis(seed)
        self.preAPFParm = self.N.getParmVal()
        self.preAPFf = self.ef()
        self.mrf.efun()
        self.postAPFParm = self.N.getParmVal()
        self.postAPFf = self.ef()
        self.APFpValue = self.get_pValue(self.postTruef, self.postAPFf, self.trueParm.size())
        self.mle = self.N.getParmVal()
        self.ml = self.N.likelihood()
        if run == 4:
          return

        self.H = numpy.matrix(self.Hessian())
        svdList = svd.svd(numpy.array(self.H))[1]
        self.precision = 0.0
        for sl in svdList:
            sl_positive = max(sl,1e-14)
            self.precision += math.log(sl_positive)
        self.N.setParm(self.saveParm)
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

    def save(self):
        return saveTCR(self)

class saveTCR(object):
    def __init__(self,tcr):
        self.run = tcr.lastrun
        self.nchan = tcr.lastnchan
        self.seed = tcr.lastseed
        self.n_trajec = tcr.lastn_trajec
        if run > 0:
            # Nuisance Fit At True
            self.preTrueParm = numpy.matrix(tcr.preTrueParm)
            self.preTruef = tcr.preTruef
            self.postTrueParm = numpy.matrix(tcr.postTrueParm)
            self.postTruef = tcr.postTruef
            self.otle = numpy.matrix(tcr.otle)
            self.otml = tcr.otml
        if run > 1:
            # Square Norm Fit
            self.preSNFParm = numpy.matrix(tcr.preSNFParm)
            self.preSNFf = tcr.preSNFf
            self.postSNFParm = numpy.matrix(tcr.postSNFParm)
            self.postSNFf = tcr.postSNFf
        if run > 2:
            # Main Parm Fit:
            global MPF
            if MPF:
                self.preMPFParm = numpy.matrix(tcr.preMPFParm)
                self.preMPFf = tcr.preMPFf
                self.postMPFParm = numpy.matrix(tcr.postMPFParm)
                self.postMPFf = tcr.postMPFf
                self.MPFpValue = tcr.MPFpValue
        if run > 3:
            # All Parm Fit
            self.preAPFParm = numpy.matrix(tcr.preAPFParm)
            self.preAPFf = tcr.preAPFf
            self.postAPFParm = numpy.matrix(tcr.postAPFParm)
            self.postAPFf = tcr.postAPFf
            self.APFpValue = tcr.APFpValue
            self.pValue = tcr.APFpValue
            self.mle = numpy.matrix(tcr.mle)
            self.ml = tcr.ml

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
        self.trueParm = self.N.getParm()
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
            vec = h.ch3ssdata(n, seed, tvec, self.trueParm, self.N.rf.fitnesslist.o(0))
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
        self.trueParm = numpy.matrix(WH.trueParm)
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

def start(seed=1, nchannels=50, modelses="ch3_101p.ses", datagenhoc="ch3ssdatagen.hoc", run=0):
    return TestCR(nchannels,seed,modelses,datagenhoc,run)

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

def onerun(seed=1, nchannels=50, n_trajectory=1):
    r = tcr # global from parrun
    cvodewrap.fs.use_fixed_step = 1.0
    tt = h.startsw()
    pc = h.ParallelContext()
    id = int(pc.id())
    print "%d start seed=%d nchannels=%d"%(id, seed, nchannels)
    run = 3
    r.compute(nchannels, seed, n_trajectory, run)
    tt = h.startsw() - tt
    print "%d finish walltime=%g seed=%d nchannels=%d n_trajectory=%d"%(id, tt, seed, nchannels, n_trajectory)
    #return value suitable for bulletin board
    result = [[tt, seed, nchannels, n_trajectory, run]]
    print 'run ', run
    if run >= 1: result += [[numpy.array(r.postTrueParm), r.postTruef]]
    if run >= 2: result += [[numpy.array(r.postSNFParm), r.postSNFf]]
    if run >= 3: result += [[numpy.array(r.postMPFParm), r.postMPFf]]
    if run >= 4: result += [[numpy.array(r.postAPFParm), r.postAPFf]]
    return result


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

def run(nruns=1,nchannels=50,modelses="ch3_101p.ses",datagenhoc="ch3ssdatagen.hoc",run=4):
    TCRs = []
    ncovers = 0
    for i in range(nruns):
        TCRi = TestCR(nchannels,i+1,modelses,datagenhoc,run)
        TCRs.append(TCRi)
        # return TCRs
        if run > 0:
            if TCRi.covers:
                print i, 'IN:', TCRi.pValue
                ncovers += 1
            else:
                print i, 'OUT:', TCRi.pValue
    if run > 0:
        print ncovers, 'covers out of', nruns
    return TCRs

def mk_tcr(modelses="ch3_101p.ses", datagen=None):
    h.load_file(modelses)
    mrf = MulRunFitHandle()
    fitfun = mrf.p.pf.generatorlist.o(0).gen.fitnesslist.o(0)
    # requirement is that datagen(nchannel, seed1, seed2, xmd) returns
    # a h.Vector of measurement values corresponding to the xmd (h.Vector)
    # time values. Note that seed2 will range from 0 to n_trajectory-1.
    # seed1 and nchannel will be passed via onerun above.
    if datagen == None:
        h.load_file("exper_data.hoc")
        datagen = h.experimentalDataGenerator
    r = repro.Repro(fitfun, datagen, 1, 1)
    global tcr
    tcr = TestCR(mrf, NrnBFiltHandle(mrf), r)
    return tcr

def parrun(nruns=0,nchannels=50, n_trajectory=1, modelses="ch3_101p.ses", datagen=None):
    mk_tcr(nchannels, modelsed, datagen)

    pc = h.ParallelContext()
    pc.runworker()
    if nruns == 0:
        nruns = int(pc.nhost())
    for i in range(nruns):
        pc.submit(onerun, i+1, nchannels, n_trajectory)
    f = open('results'+str(nchannels)+'.'+str(n_trajectory)+'.'+str(nruns), "w")
    pickle.dump(nruns, f); f.flush()
    i = 0
    while (pc.working() != 0.0):
        r = pc.pyret()
        pickle.dump(r, f); f.flush()
        i += 1
        print '%d seed %d'%(i,r[0][1])
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
