import random
import math

class NoiseGen(object):
    def __init__(self):
        self.seed = 0
        self.R = random.Random()
        self.wienerTimes = [0.0]
        self.wienerValues = [0.0]
        self.injectionTimes = [0.0]
        self.injectionValues = [0.0]
        self.TOL = 1e-7

    def TOLequal(self, a, b):
        return ((a < b + self.TOL) and (a + self.TOL > b))

    def index2old(self, times):
        # Creates an index, the same length of times
        # index = [0, None, 2, None] means:
        #        times[0] = wienerTimes[0]
        #        times[1] and times[3] are new ie not in wienerTimes
        #        times[2] = wienerTimes[2] thus wienerTimes[1] not in times
        index = []
        wi = 0
        for t in times:
            while wi < len(self.wienerTimes) and self.wienerTimes[wi] < t:
                wi += 1
            if wi < len(self.wienerTimes) and self.TOLequal(t, self.wienerTimes[wi]):
                index.append(wi)
            else:
                index.append(None)
        assert len(times) == len(index)
        return index

    def addBetween(self, unionTimes, unionValues, k, t):
        assert self.wienerTimes[k-1] < t
        assert t < self.wienerTimes[k]
        assert unionTimes[-1] < t
        if unionTimes[-1] > self.wienerTimes[k-1]:
            t0 = unionTimes[-1]
            W0 = unionValues[-1]
        else:
            t0 = self.wienerTimes[k-1]
            W0 = self.wienerValues[k-1]
        if self.TOLequal(unionTimes[-1],self.wienerValues[k-1]):
            assert self.TOLequal(unionValues[-1],self.wienerValues[k-1])
        dt0 = t - t0
        dt1 = self.wienerTimes[k] - t
        dW = self.wienerValues[k] - W0
        dt = dt0 + dt1
        mu = dW*dt0/dt
        sig = math.sqrt(dt0*dt1/dt)
        dW0 = self.R.normalvariate(mu,sig)
        unionTimes.append(t)
        unionValues.append(W0 + dW0)

    def addEnd(self, unionTimes, unionValues, t):
        assert(t > unionTimes[-1])
        assert(t > self.wienerTimes[-1])
        if unionTimes[-1] > self.wienerTimes[-1]:
            t0 = unionTimes[-1]
            W0 = unionValues[-1]
        else:
            t0 = self.wienerTimes[-1]
            W0 = self.wienerValues[-1]
        if self.TOLequal(unionTimes[-1], self.wienerTimes[-1]):
            assert self.TOLequal(unionValues[-1],self.wienerValues[-1])
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
                    unionTimes.append(self.wienerTimes[k])
                    unionValues.append(self.wienerValues[k])
                    if k == index[i]:
                        currentValues.append(unionValues[-1])
                        # print 'CV already', len(currentValues)
                    k += 1
            else:
                while k < len(self.wienerTimes) and self.wienerTimes[k] < times[i]:
                    unionTimes.append(self.wienerTimes[k])
                    unionValues.append(self.wienerValues[k])
                    k += 1
                if k >= len(self.wienerTimes):
                    self.addEnd(unionTimes, unionValues, times[i])
                    # print 'CV end', len(currentValues) + 1
                else:
                    self.addBetween(unionTimes, unionValues, k, times[i])
                    # print 'CV between', len(currentValues) + 1
                currentValues.append(unionValues[-1])
        while k < len(self.wienerTimes):
            unionTimes.append(self.wienerTimes[k])
            unionValues.append(self.wienerValues[k])
            k += 1
        return (unionTimes, unionValues, currentValues)

    def refine(self, times):
        times.sort()
        assert self.TOLequal(times[0],0.0)
        index = self.index2old(times)
        assert len(index) == len(times)
        (self.wienerTimes, self.wienerValues, self.injectionValues) = self.constructUnion(times,index)
        self.injectionTimes = times
        assert len(self.injectionTimes) == len(self.injectionValues)

    def getTW(self):
        return (self.wienerTimes, self.wienerValues)

    def getTI(self):
        return (self.injectionTimes, self.injectionValues)