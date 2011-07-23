import scipy.linalg
import numpy
import random

def ch3HMA(q1=.5, q2=0.25, nchannels=1):
    smallQ = numpy.matrix([[-q2, q2, 0],[0, -q1, q1],[0, 0, 0]])
    pX0 = numpy.matrix([[1, 0, 0]])
    states_c = numpy.matrix([[1],[1],[0]])
    A = HMA(smallQ,pX0,states_c,nchannels)
    return A

class HMA(object):
    def __init__(self,smallQ,pX0,states_c,nchannels=1):
        self.smallQ = smallQ
        self.pX0 = pX0
        self.states_c = states_c
        self.nchannels = nchannels
    
    #def sim(self,seeds,Delta,N):
     
class HMA3(object):
    def __init__(self, q1=.5, q2=.25, nchannels=1):
        self.q2 = q2
        self.q1 = q1
        self.smallQ = numpy.matrix([[-q2, q2, 0],[0, -q1, q1],[0, 0, 0]])
        self.pX0 = numpy.matrix([[1,0,0]])
        self.states_c = numpy.matrix([[1],[1],[0]])
        self.nchannels = nchannels
        self.R = random.Random()
        
    def sim(self, seeds=0, Delta=0.1, tstop=20):
        self.seeds = seeds
        self.Delta = Delta
        self.tstop = tstop
        self.N = numpy.floor(tstop/Delta)
        self.R.seed(seeds)
        obs = []
        for i in range(self.nchannels):
            t2 = self.R.expovariate(self.q2)
            t1 = self.R.expovariate(self.q1)
            ob = int(numpy.ceil((t2+t1)/Delta))
            if ob > self.N:
                ob = self.N
            obs.append(ob)
        self.simData = obs
    
    def like(self, dataObject):
        
        