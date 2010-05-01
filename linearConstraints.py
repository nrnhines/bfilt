import numpy
from numpy import matlib as matlib

class LinearConstraints(object):
    def __init__(self, nstates):
        self.equal_Dd= []
        self.greater_Dd = []
        self.nstates = nstates
        self.great0 = matlib.zeros([1,nstates])
        self.less1 = matlib.zeros([1,nstates])

    def pickStates(self,states):
        Di = matlib.zeros([1,self.nstates])
        for j in states:
            Di[0,j] = 1.0
        return Di

    def sumStatesEqual1(self,states):
        assert(len(states) > 1)
        Di = self.pickStates(states)
        di = 1.0

    def stateGreaterThan0(self,state):
        assert(self.great0[0,state] == 0)
        Di = self.pickStates([state])
        di = 0.0
        self.greater_Dd.append([Di,di])
        self.gre

    def stateLessThan1(self,state):
        assert(self.less1[0,state] == 0)
        Di = -self.pickStates([state])
        di = -1.0
        self.greater_Dd.append([Di,di])