import noise
import numpy
import random
import math
from myscipy import linalg
from neuron import h
import fitglobals
import eve
import copy
import cvodewrap

class StochasticModel(object):
    def __init__(self, dim, tlast):
        self.hhB = False
        self.nNa = 1000
        self.nK = 1000
        self.states = h.Vector()
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
        if self.hhB:
            cvodewrap.states(self.states)
            h.rates_hh(self.states[0])
            alpha_m = h.mtau_hh * h.minf_hh
            beta_m = h.mtau_hh * (1.0 - h.minf_hh)
            alpha_h = h.htau_hh * h.hinf_hh
            beta_h = h.htau_hh * (1.0 - h.hinf_hh)
            alpha_n = h.ntau_hh * h.ninf_hh
            beta_n = h.ntau_hh * (1.0 - h.ninf_hh)
            state_m = copy.deepcopy(self.states[1])
            state_h = copy.deepcopy(self.states[2])
            state_n = copy.deepcopy(self.states[3])
            #if state_m < 0.0:
            #   state_m = 0.0
            #if state_h < 0.0:
            #   state_h = 0.0
            #if state_n < 0.0:
            #    state_n = 0.0
            #if state_m > 1.0:
            #    state_m = 1.0
            #if state_h > 1.0:
            #    state_h = 1.0
            #if state_n > 1.0:
            #    state_n = 1.0
            print 'alpha_h', alpha_h, 'beta_h', beta_h, 'state_h', state_h, 'nNa', self.nNa
            self.B[0,0] = 0.0
            self.B[1,1] = math.sqrt((alpha_m*(1.0-state_m) + beta_m*state_m)/self.nNa)
            print 'Inside Sqrt', ((alpha_h*(1.0-state_h) + beta_h*state_h)/self.nNa)
            self.B[2,2] = math.sqrt((alpha_h*(1.0-state_h) + beta_h*state_h)/self.nNa)
            self.B[3,3] = math.sqrt((alpha_n*(1.0-state_n) + beta_n*state_n)/self.nK)
            return self.scale*self.B*math.sqrt(tFinal - inject0)
        else:
            assert(tFinal-inject0>=0.0)
            assert(not numpy.isnan(self.scale))
            R = self.scale*self.B*math.sqrt(tFinal - inject0)
            assert(not numpy.isnan(R).any())
            return R

