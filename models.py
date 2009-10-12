import noise
import numpy
import random

class EventTimed:
    def __init__(self):
        self.Times = None
    def erange(self, tstart, tstop, dt):
        self.Times = numpy.arange(tstart+0.0, tstop+dt/2.0,  dt+0.0).tolist();
    def set(self, times):
        self.Times = times

# Parent class to subclasses
class ObserveState0(noise.Gauss):
    def __init__(self, p, i):
        # SAME AS IN GAUSS
        self.R = random.Random()
        # save local copies of P and I
        self.P = p
        self.I = i
        # calculate noise
        self.X = self.calc()
        # NEW FROM HERE ON
        # set up event times
        self.Times = EventTimed()
        self.sigma = 0.0001
        
    def mean(self, time, state, param):
        return state[0, 0]
        
    def meas(self, time, state, param):
        return self.mean(time, state, param) + self.sigma*self.eval(time)
        
    def Dstate(self, time, state, param):
        J = numpy.matrix(numpy.zeros((1, state.shape[0])))
        J[0, 0] = 1
        return J
        
    def Dnoise(self, time, state, param):
        return numpy.matrix(self.sigma)

# This class defines a list of the observation classes employed.
# All have NoiseParams P, instance i0 (for unique random seeds)
# Finally c is a python list (e.g. [o1,o2,o3]) of the observation objects.
class ObservationModel:
    def __init__(self, p, i0, c):
        self.P = p
        self.I = i0
        self.D = len(c)
        self.C = c
        k = 0
        while k < self.D:
            self.C[k].__init__(p, i0+k)
            k += 1
            
    # Change parameters
    def change(self, p):
        self.__init__(p, self.I, self.C)
        
    # Evaluate the noise at a given time
    def eval(self, time,  List=None):
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:
            E.append(self.C[k].eval(time))
        return E

    # Evaluate the mean measurement at a given time
    def mean(self, time, state, param, List=None):
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:        
            E.append(self.C[k].mean(time, state, param).tolist())
        return numpy.matrix(E)

    # Evaluate the measurement at a given time
    def meas(self, time, state, param, List=None):
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:
            E.append(self.C[k].meas(time, state, param).tolist())
        return numpy.matrix(E)
        
    def Dstate(self, time, state, param, List=None):
        if List == None:
            List = range(0, self.D)
        E = numpy.matrix(numpy.zeros((len(List), state.shape[0])))
        Eindex = 0
        for k in List:
            E[Eindex, :] = self.C[k].Dstate(time, state, param)
            Eindex += 1
        return E
        
    def Dnoise(self, time, state, param, List=None):
        if List == None:
            List = range(0, self.D)
        E = numpy.matrix(numpy.zeros((len(List), len(List))))
        Eindex = 0
        for k in List:
            E[Eindex, Eindex] = self.C[k].Dnoise(time, state, param)
            Eindex += 1
        return numpy.matrix(E) 
        

class DecayModel:
    def __init__(self, p, i0, d):
        self.P = p
        self.I = i0
        self.D = d
        self.C = GaussVector(p, i0, d)
        
    # Change parameters
    def change(self, p):
        self.__init__(p, self.I, self.D)
        
    # Evaluate the noise at a given time
    def eval(self, time):
        E = []
        for k in List:
            E.append(self.C[k].eval(time))
        return E

    # Evaluate the mean vector field at a given time
    def vfield(self, time, state,  param, discrete=None):
        return -param.A*state
    
    def flow(self, Times, state0, param, discrete=None):
        EndTime = Times[len(Times)-1]
        if len(state) == 1:
            return (math.exp(-param.A*(EndTime-Times[0])))*state0
        else:
            return (scipy.linalg.expm(-param.A*(EndTime-Times[0])))*state0
    
    def Dstate(self, Times, state0, param, discrete=None):
        EndTime = Times[len(Times)-1]
        if len(state) == 1:
            return (math.exp(-param.A*(EndTime-Times[0])))
        else:
            return (scipy.linalg.expm(-param.A*(EndTime-Times[0])))
            
    def Dnoise(self, Times, state, param, discrete=None):
 
