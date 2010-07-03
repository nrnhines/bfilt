import noise
import numpy
import random
import math
from myscipy import linalg
from neuron import h
import fitglobals
import eve

class StochasticModel(object):
    def __init__(self, dim, tlast, covGrowthTime, varTerm):
        self.B = numpy.matrix(numpy.zeros((dim,dim)))
        self.B[0,0] = 1
        self.tlast = tlast
        times = numpy.arange(0.0,tlast,1.0)
        self.Injection = eve.EventTimed(times)
        self.InitialCov = numpy.eye(dim)
        for i in range(dim):
            self.InitialCov[i,i] = Sdiag[i]**2.0

    def updateInitial(self, Sdiag):
        dim = self.B.shape[0]
        self.InitialCov = numpy.eye(dim)
        for i in range(dim):
            self.InitialCov[i,i] = Sdiag[i]**2.0

    def updateInjectionInterval(self,dt):
        self.Injection.erange(0.0,self.tlast,dt)

    def noiseJac(self, injectionTimes):
        inject0 = injectionTimes[0]
        tFinal = injectionTimes[-1]
        return self.B*math.sqrt(tFinal - inject0)
