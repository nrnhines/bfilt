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
        
    def mean(self, time, state):
        return state[0, 0]
        
    def meas(self, time, state):
        return self.mean(time, state) + self.sigma*self.eval(time)
        
    def Dstate(self, time, state):
        J = numpy.matrix(numpy.zeros((1, state.shape[0])))
        J[0, 0] = 1
        return J
        
    def Dnoise(self, time, state):
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
    def mean(self, time, state, List=None):
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:        
            E.append(self.C[k].mean(time, state).tolist())
        return numpy.matrix(E)

    # Evaluate the measurement at a given time
    def meas(self, time, state, List=None):
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:
            E.append(self.C[k].meas(time, state).tolist())
        return numpy.matrix(E)
        
    def Dstate(self, time, state, List=None):
        if List == None:
            List = range(0, self.D)
        E = numpy.matrix(numpy.zeros((len(List), state.shape[0])))
        Eindex = 0
        for k in List:
            E[Eindex, :] = self.C[k].Dstate(time, state)
            Eindex += 1
        return E
        
    def Dnoise(self, time, state, List=None):
        if List == None:
            List = range(0, self.D)
        E = numpy.matrix(numpy.zeros((len(List), len(List))))
        Eindex = 0
        for k in List:
            E[Eindex, Eindex] = self.C[k].Dnoise(time, state)
            Eindex += 1
        return numpy.matrix(E) 
        

class DecayModel:
    def __init__(self, p, i0, d):
        self.P = p
        self.I = i0
        self.D = d
        self.C = noise.GaussVector(p, i0, d)
        self.Injection = EventTimed()
        
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
    def vfield(self, time, state, discrete=None):
        return -self.P.A*state
    
    def flow(self, Times, state0, discrete=None):
        EndTime = Times[len(Times)-1]
        if len(state) == 1:
            return (math.exp(-self.P.A*(EndTime-Times[0])))*state0
        else:
            return (scipy.linalg.expm(-self.P.A*(EndTime-Times[0])))*state0
    
    def Dstate(self, Times, state0, discrete=None):
        EndTime = Times[len(Times)-1]
        if len(state) == 1:
            return (math.exp(-self.P.A*(EndTime-Times[0])))
        else:
            return (scipy.linalg.expm(-self.P.A*(EndTime-Times[0])))
            
    def Dnoise(self, Times, state0, discrete=None):
        Dn = numpy.matrix(zeros(len(state), len(Times)-1))
        if len(state) == 1:
            for Index in range(0, len(Times)-1):
                dtLast = Times[Index+1] - Times[Index]
                dt2End = Times[len(Times)-1] - Times[Index+1]
                Dn[:, Index] = math.sqrt(dtLast)*math.exp(-self.P.A*(dt2End))*self.P.B
        else:
            for index in range(0, len(Times)-1):
                dtLast = Times[Index+1] - Times[Index]
                dt2End = Time[len(Times)-1] - Times[Index+1]
                Dn[:, Index] = math.sqrt(dtLast)*scipy.linalg.expm(-self.P.A*(dt2End))*self.P.B

class Model:
    def __init__(self, sys, obs, p):
        self.Sys = sys
        self.Obs = obs
        self.P = p
        self.Sys.change(p)
        self.Obs.change(p)
        FitEvents = self.Tabulate
        
    def change(self, p):
        self.__init__(self.Sys, self.Obs, p)
        
    def Tabulate(self):
        EventSet = set(self.Sys.Injection)
        for Index in range(0, self.Obs.D):
            EventSet = EventSet | set(self.Obs.C[Index].Times)
        ObserveIndicies = [0]*self.Obs.D
        NoiseTie
