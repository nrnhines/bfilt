import fitglobals
import noise
import models
import fitEKF
import numpy
import copy

class twoDecay(object):
	def __init__(self):
		fitglobals.debugoff()
		P = noise.NoiseParams()
		P.tstop = 20
		P.dt = 0.1
		A0= 0.5
		A1= 0.1
		P.A = numpy.matrix([[A0, 0], [0,A1]])
		# P.B, next line define noise injected to each component, uncorrelated
		P.B = numpy.matrix([[0.5, 0], [0, 0.4]])
		P.InitialCov = numpy.matrix([[1,0],[0,1]])
		elist = numpy.arange(0.1,20,0.1).tolist()
		elist2 = numpy.arange(0.05,20.0,0.05).tolist()
		O1 = models.ObserveState0(P,5)
		O2 = models.ObserveStateSum(P,6)
		# O1.Times.set([2,4,6,8,10,12,14,16,18,20])
		# O2.Times.set([3,6,9,12,15,18])
		O1.Times.set(elist)
		O2.Times.set(elist)
		O1.sigma = 0.001
		O2.sigma = 0.0001
		Obs = models.ObservationModel(P,5,[O1,O2])
		Sys = models.DecayModel(P,0,2)
		# Sys.Injection.set([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
		Sys.Injection.set(elist2)
		Initial = numpy.matrix([[1.0],[2.0]])
		self.M = models.Model(Sys,Obs,P,Initial)
		self.sim(False)
		
	def sim(self, withnoise = True):
		self.Msim = copy.deepcopy(self.M)
		if withnoise:
			print "Sim with Noise left alone"
		else:
			print "Sim with Noise Turned OFF"
			# Turn noise off
			self.Msim.P.B = numpy.matrix([[0.,0.],[0.,0.]])
			self.Msim.Obs.C[0].sigma = 0
			self.Msim.Obs.C[1].sigma = 0
		self.Data = self.Msim.sim()

	def loglike(self, A0=None,A1=None):
		simParams = self.getParams()
		if A0 == None:
			A0 = simParams[0]
		if A1 == None:
			A1 = simParams[1]
		Mtest = copy.deepcopy(self.M)
		self.setParams(A0,A1,Mtest)
		LL = fitEKF.ekf(self.Data, Mtest)
		return LL
		
	def setParams(self,A0set,A1set,Mset):
		Mset.P.A = numpy.matrix([[A0set, 0], [0,A1set]])
		
	def getParams(self, Mget=None):
		if Mget == None:
			Mget = self.Msim
		A0get = Mget.P.A[0,0]
		A1get = Mget.P.A[1,1]
		return [A0get, A1get]

class simpleDecay(twoDecay):
	def __init__(self):
		P = noise.NoiseParams()
		P.tstop = 20
		P.dt = 0.1
		P.A = numpy.matrix(0.5)
		P.B = numpy.matrix([0.1])
		P.InitialCov = numpy.matrix(1)
		O1 = models.ObserveState0(P,5)
		elist = numpy.arange(0.1,20,0.1).tolist()
		elist2 = numpy.arange(0.05,20.0,0.05)
		# O1.Times.set([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
		O1.Times.set([2,4,6,8,10,12,14,16,18,20])
		# O1.Times.set(elist)
		O1.sigma = 0.001
		Obs = models.ObservationModel(P,5,[O1])
		Sys = models.DecayModel(P,0,1)
		# Sys.Injection.set(elist2)
		Sys.Injection.set([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
		self.M = models.Model(Sys,Obs,P,numpy.matrix(2))
		self.sim(False)

	def sim(self, withnoise = True):
		self.Msim = copy.deepcopy(self.M)
		if withnoise:
			print "Sim with Noise left alone"
		else:
			print "Sim with Noise Turned OFF"
			# Turn noise off
			self.Msim.P.B = numpy.matrix([0.])
			self.Msim.Obs.C[0].sigma = 0
		self.Data = self.Msim.sim()

	def setParams(self,A0set,A1set,Mset):
		Mset.P.A = numpy.matrix([A0set])

	def getParams(self, Mget=None):
		if Mget == None:
			Mget = self.Msim
		A0get = Mget.P.A[0,0]
		A1get = 0
		return [A0get, A1get]
