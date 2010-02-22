import noise
import numpy
import random
import math
from myscipy import linalg
from neuron import h
import fitglobals

toler = 1e-8

def rmDupsWithToler(L,tol):
    lenL = len(L)
    i = 0
    while i < lenL-1:
        if L[i+1] - L[i] < tol:
            L.pop(i+1)
            lenL -= 1
        else:
            i += 1

class EventTimed:
    def __init__(self, times=None):
        self.Times = times
    
    def erange(self, tstart, tstop, dt):
        self.Times = numpy.arange(tstart+0.0, tstop+dt/2.0,  dt+0.0).tolist()
    
    def set(self, times):
        self.Times = times
    
    def makeIncreasing(self, dt=0.0):
        # 1. Eliminates duplicates
        # 2. Eliminiates non-increasing times.
        # 3. Eliminates negative times
        listReturn = numpy.asarray(self.Times).tolist()
        lenList = len(listReturn)
        index = 1
        while index < lenList:
            if (listReturn[index] - listReturn[index-1]) < dt/2.0 + toler:
                del listReturn[index]
                lenList -= 1
            else:
                index += 1
        # Get rid of non-positive numbers at beginning of list
        while listReturn[0] < 0.:
            listReturn.pop(0)
        return listReturn

class EventTable:
    def __init__(self, sto, obs):
        self.Sto = sto
        self.Obs = obs
        self.tab()
    
    def newInjectionInterval(self,dt):
        print 'New II', dt
        self.Sto.updateInjectionInterval(dt)
        self.tab()
    
    def tab(self):
        # Create collectionTimes list
        self.collectionTimes = []
        for Index in range(0, self.Obs.D):
            self.collectionTimes += self.Obs.C[Index].Times.makeIncreasing()
        self.collectionTimes.sort()
        rmDupsWithToler(self.collectionTimes,toler)
        
        # Create injectionTimes list
        self.injectionTimes = []
        Inj = self.Sto.Injection.makeIncreasing()
        # Inj.sort() # makeincreasing does this already -- remove either here or above
        # rmDupsWithToler(self.collection.TImes,toler)
        injectionsForThisData = [0]
        for timeNextOb in self.collectionTimes:
            while len(Inj) > 0 and Inj[0] < timeNextOb + toler:
                injectionsForThisData.append(Inj[0])
                Inj.pop(0)
            if len(injectionsForThisData) == 1:
                injectionsForThisData.append(timeNextOb)
            self.injectionTimes.append(injectionsForThisData)
            injectionsForThisData = [injectionsForThisData[-1]]  # start next list with last item of previous
        print 'ITs', self.injectionTimes
        
        #Create ObsNum list
        self.ObsNum = []
        OTimes = []  # OTimes will be a list of lists of event times, each list for a different obs
        for ObNum in range(0, self.Obs.D):
            OTimes.append(self.Obs.C[ObNum].Times.makeIncreasing())
        for timeNextOb in self.collectionTimes:
            ObsThisData = []
            for ObNum in range(0, self.Obs.D):
                if len(OTimes[ObNum]) > 0 and OTimes[ObNum][0] < timeNextOb + toler:
                    ObsThisData.append(ObNum)
                    temp = OTimes[ObNum].pop(0)  # delete from list
            assert(len(ObsThisData)>0)
            self.ObsNum.append(ObsThisData)
