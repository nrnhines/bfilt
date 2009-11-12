from neuron import h
h.load_file('pas.ses')

import noise
import models
import fitEKF
import numpy
import fitglobals

fitglobals.debugon()

P = noise.NoiseParams()
P.tstop = 20
P.dt = 0.1
P.A = numpy.matrix(0.5)
P.B = numpy.matrix([0.5, 0.4, 0.3])
P.InitialCov = numpy.matrix(1)
O1 = models.ObserveState0(P,5)
O2 = models.ObserveState0(P,6)
O1.Times.set([2,4,6,8,10,12,14,16,18,20])
O2.Times.set([3,6,9,12,15,18])
O1.sigma = 0.001
O2.sigma = 0.0001
Obs = models.ObservationModel(P,5,[O1,O2])
#Sys = models.DecayModel(P,0,3)
Sys = models.NeuronModel(P,0,3)
Sys.Injection.set([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
M = models.Model(Sys,Obs,P)

# Turn noise off and simulate
P.B = numpy.matrix([0.,0.,0.])
  # Not needed: change(P)
M.Obs.C[0].sigma = 0
M.Obs.C[1].sigma = 0
Data = M.sim()

# Turn noise back on and calculate log-likelihood
P.B = numpy.matrix([0.5, 0.4, 0.3])
  # Not needed: change(P)
M.Obs.C[0].sigma = 0.001
M.Obs.C[1].sigma = 0.0001
LL = fitEKF.ekf(Data,M)
