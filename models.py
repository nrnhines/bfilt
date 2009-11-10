import noise
import numpy
import random
import math
from myscipy import linalg
from neuron import h
useNeuron = False
from fitglobals import debug

class EventTimed:
    def __init__(self, times=None):
        self.Times = times
    def erange(self, tstart, tstop, dt):
        self.Times = numpy.arange(tstart+0.0, tstop+dt/2.0,  dt+0.0).tolist()
    def set(self, times):
        self.Times = times
    # The following function does several things:
    # 1. Rounds event times to a particular discretization (dt)
    # 2. After rounding, eliminates duplicates
    # 3. Eliminiates non-increasing times.
    def round(self, dt):
        matTimes = numpy.asarray(self.Times)
        roundTimes = (matTimes/(dt+0.0)).round()
        roundTimes = dt*roundTimes
        listReturn = roundTimes.tolist()
        listReturn = listReturn
        lenList = len(listReturn)
        index = 1
        while index < lenList:
            if (listReturn[index] - listReturn[index-1]) < dt/2.0:
                del listReturn[index]
                lenList -= 1
            else:
                index += 1
        # Get rid of non-positive numbers at beginning of list
        while listReturn[0] <= 0.:
            listReturn.pop(0)
        return listReturn

# Parent class to subclasses
class ObserveState0(noise.Gauss):
    def __init__(self, p, i, times=None):
        self.__init(p, i,  times)

    def __init(self, p, i, times=None):
        self.Times = EventTimed(times)
        self.sigma = 0.001
        self.observed = 0
        self.base_init(p, i)
    
    # just change P not I
    def change(self, p):
        self.base_init(p, self.I)
    
    def mean(self, time, state):
        return state[self.observed, 0]
        
    def meas(self, time, state):
        return self.mean(time, state) + self.sigma*self.eval(time)
        
    def Dstate(self, time, state):
        J = numpy.matrix(numpy.zeros((1, state.shape[0])))
        J[0, self.observed] = 1
        return J
        
    def Dnoise(self, time, state):
        return numpy.matrix(self.sigma)

class ObserveStateSum(ObserveState0):
    def __init__(self, p, i, times=None):
        self.__init(p, i,  times)

    def __init(self, p, i, times=None):
        self.Times = EventTimed(times)
        self.sigma = 0.001
        self.base_init(p, i)

    def mean(self, time, state):
        ob = 0
        for observed in range(state.shape[0]):
           ob += state[observed, 0]
        return ob
    
    def Dstate(self, time, state):
        J = numpy.matrix(numpy.zeros((1, state.shape[0])))
        for observed in range(state.shape[0]):
            J[0, observed] = 1
        return J

class NeuronObservable(ObserveState0):
    
    def __init__(self, hpointer, p, i, times=None):
       if (times == None):
         raise RuntimeError, "unrecoverable"
       self.hpt = hpointer
       ObserveState0.__init__(self, p, i, times)

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
                df = self.mean(time,state)
                x[i] = temp
                DFx[:,i] = (df - value)/h
        return DFx
   
    
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
            self.C[k].base_init(p, i0+k)
            k += 1
            
    # Change parameters
    def change(self, p):
        self.__init__(p, self.I, self.C)
        
    # Evaluate the noise at a given time
    def eval(self, time,  List=None):
        time = time[len(time)-1]
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:
            E.append(self.C[k].eval(time))
        return E

    # Evaluate the mean measurement at a given time
    def mean(self, time, state, List=None):
        time = time[len(time)-1]
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:        
            E.append(self.C[k].mean(time, state))
        return numpy.matrix(E).T

    # Evaluate the measurement at a given time
    def meas(self, time, state, List=None):
        time = time[len(time)-1]
        E = []
        if List == None:
            List = range(0, self.D)
        for k in List:
            E.append(self.C[k].meas(time, state))
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
        
class NeuronModel(object):
    def __init__(self, p, i0, d, times=None):
        self.P = p
        self.I = i0
        self.D = d
        self.V = noise.GaussVector(p, i0, d)
        self.Injection = EventTimed(times)
        h.cvode.atol(1e-6)
        h.cvode_active(1)
        h.stdinit()
      
    # This function has not changed  
    def change(self, p):
        self.__init__(p, self.I, self.D, self.Injection.Times)
        
    # This function has not changed
    def eval(self, time):
        E = numpy.matrix(numpy.zeros((self.D, 1)))
        for k in range(0, self.D):
            E[k, :] = self.V.C[k].eval(time)
        return E
        
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

    def stochflow(self, Times, state0, discrete=None):
        if discrete:
            discrete.restore()
            self.moveto(Times[0])
        assert(h.t == Times[0])
        x = numpy.matrix(state0)
        Wm = self.eval(Times[0])
        for t in Times[1:]:
            if debug:
                print 'x', x
            h.cvode.yscatter(h.Vector(x)) #MH had x.T
            h.cvode.re_init()
            h.cvode.solve(t)
            s = h.Vector()
            h.cvode.states(s)
            x =  numpy.matrix(s).T
            W = self.eval(t)
            if debug:
                print 'B', self.P.B
                print 'dW', W - Wm
                print 'x', x
            x += (self.P.B*(W - Wm))
            Wm = W
        return x
      
    def perturbedflow(self, Times, state0, iTimes, perturb, discrete=None):
        if discrete:
            discrete.restore()
            self.moveto(Times[0])
        assert(h.t == Times[0])
        x = numpy.matrix(state0)
        if debug:
            print 'x after assert', x
        for i in range(1, len(Times)):
            t = Times[i]
            h.cvode.yscatter(h.Vector(x))
            h.cvode.re_init()
            h.cvode.solve(t)
            s = h.Vector()
            h.cvode.states(s)
            x =  numpy.matrix(s).T
        if debug:
            print 'x to perturb', x
            print 'perturb to perturb', perturb
        if i == iTimes:
            x += perturb
        return x
      
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

    def Dnoise(self, Times, state0, discrete=None):
        x = numpy.matrix(state0)
        value = self.flow(Times, x, discrete)
        ncol = (len(Times) - 1) * self.D
        DFx = numpy.matrix(numpy.zeros((len(value), ncol)))
        sqrtEps = math.sqrt(numpy.finfo(numpy.double).eps)
        i = 0
        sqrtEps = 1e-3
        h = sqrtEps
        e = numpy.matrix(numpy.zeros((self.D,1)))
        for itimes in range(1, len(Times)):
            dW = math.sqrt(Times[itimes] - Times[itimes-1])
            hdW = h*dW
            for idW in range(self.D):
                e[idW, 0] = 1.0
                perturb = self.P.B*e*hdW
        #OK for now but we need to reprogram for efficiency!
                df = self.perturbedflow(Times, x, itimes, perturb, discrete)
                DFx[:,i] = (df - value)/h
                e[idW, 0] = 0.0
                i += 1
        return DFx

class DecayModel:
    def __init__(self, p, i0, d, times=None):
        self.P = p
        self.I = i0
        self.D = d
        self.V = noise.WienerVector(p, i0, d)
        self.Injection = EventTimed(times)
        
    # Change parameters
    def change(self, p):
        self.__init__(p, self.I, self.D, self.Injection.Times)
        
    # Evaluate the noise at a given time
    def eval(self, time):
        E = numpy.matrix(numpy.zeros((self.D, 1)))
        for k in range(0, self.D):
            E[k, :] = self.V.C[k].eval(time)
        return E
        
    def dim(self):
        return self.P.A.shape[1]
        
    # Evaluate the mean vector field at a given time
    def vfield(self, time, state, discrete=None):
        return -self.P.A*state
    
    def flow(self, Times, state0, discrete=None):
        EndTime = Times[len(Times)-1]
        if len(state0) == 1:
            return (math.exp(-self.P.A*(EndTime-Times[0])))*state0
        else:
            return (linalg.expm(-self.P.A*(EndTime-Times[0])))*state0

    def stochflow(self, Times, state0, discrete=None):
        for interval in range(1, len(Times)):
            state0 = self.flow([Times[interval-1], Times[interval]], state0, discrete)
            temp = state0
            state0 += self.P.B*(self.eval(Times[interval])-self.eval(Times[interval-1]))
            # print 'diff', state0-temp
            if debug:
                print 'state0', state0
                print 'final state:', state0, '@ Times =', Times
        return state0
        
    def Dstate(self, Times, state0, discrete=None):
        EndTime = Times[len(Times)-1]
        if len(state0) == 1:
            return numpy.matrix(math.exp(-self.P.A*(EndTime-Times[0])))
        else:
            return numpy.matrix(linalg.expm(-self.P.A*(EndTime-Times[0])))
            
    def Dnoise(self, Times, state0, discrete=None):
        NumNoise = self.D
        Dn = numpy.matrix(numpy.zeros((len(state0), (len(Times)-1)*NumNoise)))  
        if len(state0) == 1:
            for IndexNum in range(0, (len(Times)-1)):  # loop over left-endpoints
                dtLast = Times[IndexNum+1] - Times[IndexNum]
                dt2End = Times[len(Times)-1] - Times[IndexNum+1]
                Dn[:, IndexNum*NumNoise:(IndexNum+1)*NumNoise] = numpy.matrix(math.sqrt(dtLast)*math.exp(-self.P.A*(dt2End)))*self.P.B
        else:   # only difference use expm instead of exp
            for IndexNum in range(0, len(Times)-1):
                dtLast = Times[IndexNum+1] - Times[IndexNum]
                dt2End = Times[len(Times)-1] - Times[IndexNum+1]
                Dn[:, IndexNum*NumNoise:(IndexNum+1)*NumNoise] = numpy.matrix(math.sqrt(dtLast)*linalg.expm(-self.P.A*(dt2End)))*self.P.B
        return Dn


class Model:
    def __init__(self, sys, obs, p,  initial=None):
        self.Sys = sys
        #  if initial == None:
        #     initial = numpy.ones((self.Sys.dim(), 1), float)
        
        # TAKE THESE COMMENTS OUT FOR NEURON MODELS!!!
        # h.stdinit()
        # s = h.Vector()
        # h.cvode.states(s)
        # initial = numpy.matrix(s).T
        # print initial
        
        self.Initial = initial
        self.Obs = obs
        self.P = p
        self.Sys.change(p)
        self.Obs.change(p)
        self.FitEvents = self.Tabulate()
        
    def change(self, p):
        self.__init__(self.Sys, self.Obs, p, self.Initial)
        
    def Tabulate(self):
        Inj = self.Sys.Injection.round(self.P.dt)
        Os = [];
        for Index in range(0, self.Obs.D):
            Os.append(self.Obs.C[Index].Times.round(self.P.dt)) # append Obs times
            Inj += Os[len(Os)-1]  # concatenate last set of times
        Inj.sort()
        InjTimes = EventTimed(Inj)
        Inj = InjTimes.round(self.P.dt)
        ObsEvents = []
        InjectionsForThisData = [0]
        table = []
        for i in range(0, len(Inj)):
            # ObsEvents.append([])
            InjectionsForThisData.append(Inj[i])
            ObsEvents = []
            for ObNum in range(0, self.Obs.D):
                if len(Os[ObNum])>0 and Os[ObNum][0] < Inj[i] + (self.P.dt/2.0):
                    # ObsEvents[i].append(ObNum)
                    ObsEvents.append(ObNum)
                    temp = Os[ObNum].pop(0)
            if len(ObsEvents) > 0:
                table.append([InjectionsForThisData, ObsEvents])
                InjectionsForThisData = [Inj[i]]
        return table
        
    def sim(self):
        Data = []
        state = self.Initial
        for FEindex in range(0, len(self.FitEvents)):
            Times = self.FitEvents[FEindex][0]
            state = self.Sys.stochflow(Times, state)
            ObsNum = self.FitEvents[FEindex][1]
            Data.append(self.Obs.meas(Times, state, ObsNum))
        return Data
