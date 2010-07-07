import noise
import numpy
import random
import math
from myscipy import linalg
from neuron import h
import fitglobals
import eve

class StochasticModel(object):
    def __init__(self, dim, tlast):
        self.B = numpy.matrix(numpy.zeros((dim,dim)))
	self.scale = 1.0
        self.B[0,0] = 1
        self.tlast = tlast
        times = numpy.arange(0.0,tlast,1.0)
        self.Injection = eve.EventTimed(times)
        self.InitialCovSqrt = numpy.eye(dim)
        self.InitialCov = self.InitialCovSqrt*self.InitialCovSqrt.T

    def updateInitial(self):
        self.InitialCov = self.InitialCovSqrt*self.InitialCovSqrt.T

    def updateInjectionInterval(self,dt):
        self.Injection.erange(0.0,self.tlast,dt)

    def noiseJac(self, injectionTimes):
        inject0 = injectionTimes[0]
        tFinal = injectionTimes[-1]
        return self.scale*self.B*math.sqrt(tFinal - inject0)
