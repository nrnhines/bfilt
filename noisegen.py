import random
import math

class Wiener(object):
    def __init__(self):
        self.seed = 0
        self.R = random.Random()
        self.R.seed(self.seed)
        self.processTimes = [0.0]
        self.processValues = [0.0]
        self.evalTimes = [0.0]
        self.evalValues = [0.0]
        self.TOL = 1e-7

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

    def refine(self, times):
        times.sort()
        assert self.propertimes(times)
        index = self.index2old(times)
        assert len(index) == len(times)
        (self.processTimes, self.processValues, self.evalValues) = self.constructUnion(times,index)
        self.evalTimes = times
        assert len(self.evalTimes) == len(self.evalValues)

    def getProc(self):
        return (self.processTimes, self.processValues)

    def getEval(self):
        return (self.evalTimes, self.evalValues)

class Gauss(Wiener):
    def __init__(self):
        self.seed = 0
        self.R = random.Random()
        self.R.seed(self.seed)
        self.processTimes = []
        self.processValues = []
        self.evalTimes = []
        self.evalValues = []
        self.TOL = 1e-7

    def propertimes(self,times):
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
