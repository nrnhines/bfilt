import fitglobals
import noise
import models
import fitEKF
import numpy


fitglobals.debugon()

P = noise.NoiseParams()
P.tstop = 20
P.dt = 0.1
A0= 0.5
A1= 0.1
P.A = numpy.matrix([[A0, 0], [0,A1]])
# P.B, next line define noise injected to each component, uncorrelated
P.B = numpy.matrix([[0.5, 0], [0, 0.4]])
P.InitialCov = numpy.matrix([[1,0],[0,1]])
O1 = models.ObserveState0(P,5)
O2 = models.ObserveStateSum(P,6)
O1.Times.set([2,4,6,8,10,12,14,16,18,20])
O2.Times.set([3,6,9,12,15,18])
O1.sigma = 0.001
O2.sigma = 0.0001
Obs = models.ObservationModel(P,5,[O1,O2])
Sys = models.DecayModel(P,0,2)
Sys.Injection.set([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
Initial = numpy.matrix([[1.0],[2.0]])
M = models.Model(Sys,Obs,P,Initial)

# Turn noise off and simulate
P.B = numpy.matrix([[0.,0.],[0.,0.]])
# Not needed: change(P)
M.Obs.C[0].sigma = 0
M.Obs.C[1].sigma = 0
Data = M.sim()

# Turn noise back on and calculate log-likelihood
P.B = numpy.matrix([[0.5, 0], [0, 0.4]])
# Not needed: change(P)
M.Obs.C[0].sigma = 0.001
M.Obs.C[1].sigma = 0.0001
print 'debug', fitglobals.debug
LL = fitEKF.ekf(Data,M)

fitglobals.debugoff()

def setParams(Mset,A0set,A1set):
    Mset.P.A = numpy.matrix([[A0set, 0], [0,A1set]])

def getParams(Mget):
    A0get = Mget.P.A[0,0]
    A1get = Mget.P.A[1,1]
    return [A0get, A1get]
