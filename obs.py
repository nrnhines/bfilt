import noise
import numpy
import random
import math
from myscipy import linalg
from neuron import h
import fitglobals
import eve

# Parent class to subclasses
class ObserveStateK(object):
    def __init__(self, times=None):
        self.__init(times)
    
    def __init(self, times=None):
        self.Times = eve.EventTimed(times)
        self.sigma = 0.001
        self.K = 0
    
    def mean(self, time, state):
        return state[self.K, 0]
    
    def Dstate(self, time, state):
        J = numpy.matrix(numpy.zeros((1, state.shape[0])))
        J[0, self.K] = 1
        return J
    
    def Dnoise(self, time, state):
        return numpy.matrix(self.sigma)

class NeuronObservable(ObserveStateK):
    #TAKE OUT p and i: def __init__(self, hpointer, p, i, times=None):
    def __init__(self, hpointer, times=None):
        if (times == None):
            raise RuntimeError, "unrecoverable"
        self.hpt = hpointer
        ObserveStateK.__init__(self, times)
    
    def mean(self, time, state):  # the observable (under zero noise, ie mean)
        ss = h.Vector() #save
        h.cvode.states(ss)
        
        ds = h.Vector()
        h.cvode.f(time, h.Vector(state), ds)
        #measurement goes here
        x = self.hpt.val
        
        h.cvode.yscatter(ss) #restore
        return x
    
    def Dstate(self,time,state):  # Derivative of Observable w.r.t. final state
        x = numpy.matrix(state)
        value = self.mean(time,state)
        DFx = numpy.matrix(numpy.zeros((1, len(x))))
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
            df = self.mean(time,x)
            x[i] = temp
            DFx[:,i] = (df - value)/h
        return DFx


# This class defines a list of the observation classes employed.
# All have NoiseParams P, instance i0 (for unique random seeds)
# Finally c is a python list (e.g. [o1,o2,o3]) of the observation objects.
class ObservationModel:
    def __init__(self, c):
        self.D = len(c)
        self.C = c
    
    # Evaluate the mean measurement at a given time
    def mean(self, time, state, List=None):
        time = time[len(time)-1]
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:
            E.append(self.C[k].mean(time, state))
        return numpy.matrix(E).T
    
    def Dstate(self, time, state, List=None):
        time = time[len(time)-1]
        if List == None:
            List = range(0, self.D)
        E = numpy.matrix(numpy.zeros((len(List), state.shape[0])))
        Eindex = 0
        for k in List:
            E[Eindex, :] = self.C[k].Dstate(time, state)
            Eindex += 1
        return E
    
    def Dnoise(self, time, state, List=None):
        time = time[len(time)-1]
        if List == None:
            List = range(0, self.D)
        E = numpy.matrix(numpy.zeros((len(List), len(List))))
        Eindex = 0
        for k in List:
            E[Eindex, Eindex] = self.C[k].Dnoise(time, state)
            Eindex += 1
        return numpy.matrix(E)
