import random
import math
import copy
import numpy
import scipy
import scipy.linalg as linalg
import eve
import sto
import detsys
from neuron import h

class Wiener(object):
    def __init__(self, seed):
        self.seed = seed
        self.R = random.Random()
        self.R.seed(self.seed)
        (self.processTimes, self.processValues) = self.initProcess()
        (self.evalTimes, self.evalValues) = self.initProcess()
        self.TOL = 1e-7
        self.dWList = [0.0]

    def initProcess(self):
        return ([0.0],[0.0])

    def TOLequal(self, a, b):
        return ((a <= b + self.TOL) and (a + self.TOL >= b))

    def index2old(self, times):
        # Creates an index, the same length of times
        # index = [0, None, 2, None] means:
        #        times[0] = processTimes[0]
        #        times[1] and times[3] are new ie not in processTimes
        #        times[2] = processTimes[2] thus processTimes[1] not in times
        index = []
        wi = 0
        for t in times:
            while wi < len(self.processTimes) and self.processTimes[wi] < t:
                wi += 1
            if wi < len(self.processTimes) and self.TOLequal(t, self.processTimes[wi]):
                index.append(wi)
            else:
                index.append(None)
        assert len(times) == len(index)
        return index

    def addBetween(self, unionTimes, unionValues, k, t):
        assert self.processTimes[k-1] < t
        assert t < self.processTimes[k]
        assert unionTimes[-1] < t
        if unionTimes[-1] > self.processTimes[k-1]:
            t0 = unionTimes[-1]
            W0 = unionValues[-1]
        else:
            t0 = self.processTimes[k-1]
            W0 = self.processValues[k-1]
        if self.TOLequal(unionTimes[-1],self.processTimes[k-1]):
            assert self.TOLequal(unionValues[-1],self.processValues[k-1])
        dt0 = t - t0
        dt1 = self.processTimes[k] - t
        dW = self.processValues[k] - W0
        dt = dt0 + dt1
        mu = dW*dt0/dt
        sig = math.sqrt(dt0*dt1/dt)
        dW0 = self.R.normalvariate(mu,sig)
        unionTimes.append(t)
        unionValues.append(W0 + dW0)

    def addEnd(self, unionTimes, unionValues, t):
        assert(t > unionTimes[-1])
        assert(t > self.processTimes[-1])
        if unionTimes[-1] > self.processTimes[-1]:
            t0 = unionTimes[-1]
            W0 = unionValues[-1]
        else:
            t0 = self.processTimes[-1]
            W0 = self.processValues[-1]
        if self.TOLequal(unionTimes[-1], self.processTimes[-1]):
            assert self.TOLequal(unionValues[-1],self.processValues[-1])
        dt = t - t0
        mu = 0.0
        sig = math.sqrt(dt)
        dW0 = self.R.normalvariate(mu,sig)
        unionTimes.append(t)
        unionValues.append(W0 + dW0)

    def constructUnion(self, times, index):
        k = 0
        unionTimes = []
        unionValues = []
        currentTimes = []
        currentValues = []
        for i in range(len(index)):
            # print 'k', k, 'i', i, 'index[i]', index[i]
            if not index[i] == None:
                while k <= index[i]:
                    unionTimes.append(self.processTimes[k])
                    unionValues.append(self.processValues[k])
                    if k == index[i]:
                        currentValues.append(unionValues[-1])
                        # print 'CV already', len(currentValues)
                    k += 1
            else:
                while k < len(self.processTimes) and self.processTimes[k] < times[i]:
                    unionTimes.append(self.processTimes[k])
                    unionValues.append(self.processValues[k])
                    k += 1
                if k >= len(self.processTimes):
                    self.addEnd(unionTimes, unionValues, times[i])
                    # print 'CV end', len(currentValues) + 1
                else:
                    self.addBetween(unionTimes, unionValues, k, times[i])
                    # print 'CV between', len(currentValues) + 1
                currentValues.append(unionValues[-1])
        while k < len(self.processTimes):
            unionTimes.append(self.processTimes[k])
            unionValues.append(self.processValues[k])
            k += 1
        return (unionTimes, unionValues, currentValues)

    def propertimes(self,times):
        return self.TOLequal(times[0],0.0)

    def dWeval(self):
        self.dWList = [0.0]
        for k in range(1,len(self.evalValues)):
            self.dWList.append(self.evalValues[k] - self.evalValues[k-1])

    def refine(self, times):
        times.sort()
        assert self.propertimes(times)
        index = self.index2old(times)
        assert len(index) == len(times)
        (self.processTimes, self.processValues, self.evalValues) = self.constructUnion(times,index)
        self.evalTimes = times
        assert len(self.evalTimes) == len(self.evalValues)
        self.dWeval()

    def reseed(self,seed):
        self.seed = seed
        self.R.seed(seed)
        (self.processTimes, self.processValues) = self.initProcess()
        self.refine(self.evalTimes)

    def getProc(self):
        return (self.processTimes, self.processValues)

    def getEval(self):
        return (self.evalTimes, self.evalValues)

class Gauss(Wiener):
    def __init__(self):
        self.seed = 10000
        self.R = random.Random()
        self.R.seed(self.seed)
        (self.processTimes, self.processValues) = self.initProcess()
        (self.evalTimes, self.evalValues) = self.initProcess()
        self.dWList = None
        self.TOL = 1e-7

    def initProcess(self):
        return ([], [])

    def dWeval(self):
        self.dWList = None

    def propertimes(self,times):
        if len(times) == 0:
            return True
        else:
            return times[0] >= 0.0

    def addBetween(self, unionTimes, unionValues, k, t):
        assert self.processTimes[k-1] < t
        assert t < self.processTimes[k]
        assert unionTimes[-1] < t
        if self.TOLequal(unionTimes[-1],self.processTimes[k-1]):
            assert self.TOLequal(unionValues[-1],self.processValues[k-1])
        mu = 0
        sig = 1
        G = self.R.normalvariate(mu,sig)
        unionTimes.append(t)
        unionValues.append(G)

    def addEnd(self, unionTimes, unionValues, t):
        bigUnion = (len(unionTimes) > 0)
        bigProc = (len(self.processTimes) > 0)
        if bigUnion:
            assert(t > unionTimes[-1])
        if bigProc:
            assert(t > self.processTimes[-1])
        if bigUnion and bigProc and self.TOLequal(unionTimes[-1], self.processTimes[-1]):
            assert self.TOLequal(unionValues[-1],self.processValues[-1])
        mu = 0
        sig = 1
        G = self.R.normalvariate(mu,sig)
        unionTimes.append(t)
        unionValues.append(G)

class WienerVector(object):
    def __init__(self,dim,seed0):
        self.seed0 = seed0
        self.dim = dim
        self.offset = 10000
        self.C = []
        for k in range(dim):
            self.C.append(Wiener(self.seed0+k*self.offset))

    def reseed(self,seed0):
        self.seed0 = seed0
        for k in range(self.dim):
            self.C[k].reseed(self.seed0+k*self.offset)

    def refine(self,times):
        for k in range(self.dim):
            self.C[k].refine(times)

    def dW(self,k):
        dWk = []
        for i in range(self.dim):
            dWk.append(self.C[i].dWList[k])
        dWk = numpy.matrix(dWk)
        return dWk.T

class Initial(object):
    def __init__(self, n):
        self.n = n
        self.R = random.Random()
        self.seed = 20000
        self.g = []
        self.draw() # self.R.seed set in self.draw()

    def draw(self):
        mu = 0
        sig = 1
        self.R.seed(self.seed)
        self.g = []
        for i in range(self.n):
            self.g.append(self.R.normalvariate(mu,sig))

    def gmatrix(self):
        return numpy.matrix(self.g).T

    def reseed(self,seed):
        self.seed = seed
        self.draw()  # self.R.seed set in self.draw()

    def resize(self,n):
        self.n = n
        self.draw()

    def ic(self, m, cov):
        B = linalg.cholesky(cov)
        gm = self.gmatrix()
        return m + B*gm

class Gen(object):
    def __init__(self, N):
        self.N = N
        dimWV = self.N.Eve.Sto.B.shape[1]
        self.W = WienerVector(dimWV,0)
        self.G = Gauss()
        self.IC = Initial(self.N.Sys.Initial.size)
        self.seed = 0
        self.offset = 10000
        self.reseed(self.seed)
        self.TOL = 1e-7
        self.measNum = 0

    def TOLequal(self, a, b):
        return ((a <= b + self.TOL) and (a + self.TOL >= b))

    def reseed(self,seed):
        self.seed = seed
        self.W.reseed(seed)
        self.G.reseed(seed+self.offset)
        self.IC.reseed(seed+2*self.offset)

    def eventData(self):
        TOL = 1e-7
        collectionList = copy.deepcopy(self.N.Eve.collectionTimes)
        injectionList = [0.0]
        lastInjection = 0.0
        for interval in self.N.Eve.injectionTimes:
            for it in interval:
                if it > lastInjection + TOL:
                    lastInjection = it
                    injectionList.append(it)
        ni = len(injectionList)
        ki = 0
        nc = len(collectionList)
        kc = 0
        dWi = 0
        dWindex = []
        stop = []
        iscollect = []
        while (ki < ni) or (kc < nc):
            if ki == ni:
                stop.append(collectionList[kc])
                dWindex.append(None)
                iscollect.append(True)
                kc += 1
            elif kc == nc:
                stop.append(injectionList[ki])
                dWindex.append(dWi)
                dWi += 1
                iscollect.append(False)
                ki += 1
            elif self.TOLequal(collectionList[kc], injectionList[ki]):
                stop.append(injectionList[ki])
                dWindex.append(dWi)
                dWi += 1
                iscollect.append(True)
                kc += 1
                ki += 1
            elif collectionList[kc] < injectionList[ki]:
                stop.append(collectionList[kc])
                dWindex.append(None)
                iscollect.append(True)
                kc += 1
            else:
                assert(injectionList[ki] < collectionList[kc])
                stop.append(injectionList[ki])
                dWindex.append(dWi)
                dWi += 1
                iscollect.append(False)
                ki += 1
        return (stop,dWindex,iscollect,injectionList,collectionList)

    def collect(self,time,state,k):
        m = self.N.Eve.Obs.C[self.measNum].mean(time,state)
        r = self.N.Eve.Obs.C[self.measNum].sigma*self.G.evalValues[k]
        return m + r

    def inject(self,time,state,k):
        B = self.N.Eve.Sto.B
        dW = self.W.dW(k)
        R = B*dW
        state = state + R
        return state

    def datasim(self):
        # ONEOBS ONLY
        (times,data) = self.sim()
        Data = []
        for d in data:
            Data.append(numpy.matrix(d))
        return Data

    def sim(self):
        # ONEOBS ONLY
        (stop,dWindex,iscollect,injectionList,collectionList) = self.eventData()
        self.W.refine(injectionList)
        self.G.refine(collectionList)
        state = numpy.matrix(self.IC.ic(self.N.Sys.Initial, self.N.Eve.Sto.InitialCov))
        data = []
        GIndex = 0
        if iscollect[0] and self.TOLequal(stop[0],0.0):
            data.append(self.collect(0.0,state,GIndex))
            GIndex = 1
        for k in range(1,len(stop)):
            state = numpy.matrix(self.N.Sys.flow([stop[k-1],stop[k]],state))
            if dWindex[k]:
                state = self.inject(stop[k],state,dWindex[k])
            if iscollect[k]:
                data.append(self.collect(stop[k],state,GIndex))
                GIndex += 1
        return (collectionList,data)
