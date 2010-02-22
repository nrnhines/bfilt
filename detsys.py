import noise
import numpy
import random
import math
from myscipy import linalg
from neuron import h
import fitglobals

class NeuronModel(object):
    def __init__(self):
        h.cvode.atol(1e-6)
        h.cvode_active(1)
        h.stdinit()
        s = h.Vector()
        h.cvode.states(s)
        self.Initial = numpy.matrix(s).T
    
    def dim(self):
        s = h.Vector()
        h.cvode.states(s)
        return len(s)
    
    def vfield(self, time, state, discrete=None):
        s = h.Vector(state)
        d = h.Vector()
        h.cvode.f(time, s, d)
        return numpy.matrix(d)
    
    def moveto(self, t0):
        h.t = t0
        return
        if h.t == t0:
            return
        elif t0 == 0:
            stdinit()
        elif h.t < t0:
            h.cvode.solve(t0)
        elif h.t > t0:
            h.stdinit()
            h.cvode.solve(t0)
    
    def flow(self, Times, state0, discrete=None):
        if discrete:
            discrete.restore()
        self.moveto(Times[0])
        h.cvode.yscatter(h.Vector(state0))
        h.cvode.re_init()
        assert(h.t == Times[0])
        h.cvode.solve(Times[-1])
        s = h.Vector()
        h.cvode.states(s)
        return numpy.matrix(s).T
    
    # jacobian of the flow with respect to state variables
    def Dstate(self, Times, state0, discrete=None):
        x = numpy.matrix(state0)
        value = self.flow(Times, x, discrete)
        DFx = numpy.matrix(numpy.zeros((len(value), len(x))))
        sqrtEps = math.sqrt(numpy.finfo(numpy.double).eps)
        sqrtEps = 1e-3
        for i in range(len(x)):
            temp = x[i,0]
            if abs(temp) > 1:
                h = sqrtEps*abs(temp)
            else:
                h = sqrtEps
            x[i] = temp + h
            h = x[i] - temp
            df = self.flow(Times, x, discrete)
            x[i] = temp
            DFx[:,i] = (df - value)/h
        return DFx
    
    def flowJac(self, tStart, injectionTimes,m,discrete=None):
        inject0 = injectionTimes[0]
        tFinal = injectionTimes[-1]
        mb = self.flow([tStart,tFinal],m,discrete)
        A = self.Dstate([tStart,tFinal],m,discrete)
        # MOVED TO STO: B = self.Sto.B*math.sqrt(tFinal - inject0)
        return (mb, A, tFinal)
