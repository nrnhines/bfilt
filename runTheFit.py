import noise
import models
import fitEKF
import numpy

P = noise.NoiseParams()
P.tstop = 20
P.dt = 0.1
P.A = numpy.matrix(0.5)
P.B = numpy.matrix([0.5, 0.4, 0])
P.InitialCov = numpy.matrix(1)
O1 = models.ObserveState0(P,5)
O2 = models.ObserveState0(P,6)
O1.Times.set([2,4,6,8,10,12,14,16,18,20])
O2.Times.set([3,6,9,12,15,18])
Obs = models.ObservationModel(P,5,[O1,O2])
Sys = models.DecayModel(P,0,3)
Sys.Injection.set([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
M = models.Model(Sys,Obs,P)
Data = M.sim()
LL = fitEKF.ekf(Data,M)
