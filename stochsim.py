import random
import math
import pylab
import wiener
import flow

from neuron import h
import cvodewrap

initial_Wdt = 0.001
initial_Ndt = 0.01
initial_Mdt = 0.1
initial_tMax = 5
initial_seed = 0
initial_jumpW = 0
initial_jumpM = 10000

sigm = 0.00001
sigp = 0.00001

class simData(wiener.Wiener):
    Ndt = initial_Ndt 
    Mdt = initial_Mdt
    Y = pylab.zeros(1)
    
    def __init__(self):
        self.settMax(initial_tMax)
        self.setWdt(initial_Wdt)
        self.setseed(initial_seed)
        self.setjump(initial_jumpW,initial_jumpM)
        self.sm = sigm
        self.sp = sigp
    
    def sim(self):
        tMax = self.gettMax()
        tData = 0
        tNoise = 0
        h.stdinit()
        x = h.Vector()
        Data = []
        while tData < tMax:
            noiseTimes = [] 
            tData += self.Mdt
            if tData > tMax:
                tData = tMax
            while tNoise < min(tMax,tData):
                tNoise += self.Ndt
                noiseTimes.append(tNoise)
                if tNoise > tMax:
                    tNoise = tMax
                h.continuerun(tNoise)
                cvodewrap.states(x)
                x.x[0] += self.sp*(self.evalW(tNoise) - self.evalW(tNoise-self.Ndt))
                cvodewrap.yscatter(x)
                cvodewrap.re_init()
            Data.append([noiseTimes, x.x[0] + self.sm*self.evalM(tData)])
        return Data 

def test0(a,  x0, dt,  tMax):
    t = 0
    x = x0
    Data = [[0, x0]]
    while t<=tMax:
        t += dt
        x = x*math.exp(-a*dt)
        Data.append([t, x])
    return Data
