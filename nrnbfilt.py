from neuron import h

class NrnBFilt(object):
	def __init__(self, ho):
		self.rf = ho
    ol = []
    P = noise.NoiseParam()
    P.dt = 0.1
    vl = self.rf.yvarlist
    fl = self.rf.fitnesslist
    tlast = 0
    for i in range(len(vl)):
      tl = list(fl.o(i).xdat_)
      o = models.NeuronObservable(vl.o(i), P, i, tl)
      o.sigma = 0.0001
      ol.append(o)
      if (tlast < tl[-1]):
        tlast = tl[-1]
    P.tstop = tlast
    s = h.Vector()
    cvode.states(s)
    P.B = numpy.matrix(numpy.zeros(len(s), len(vl))
    P.B[0,0] = .1
    P.InitialCov = numpy.I(len(s))
    Obs = models.ObservationModel(P, 1000, ol)
    Sys = models.NeuronModel(P, 0, len(vl))
    Sys.Injection.erange(0.0, tlast, 1.0)
    self.M = models.Model(Sys, Obs, P)
    self.Data = __data(fl,self.M.FitEvents)
    
  def __data(self,fl,FitEvents):
		# counter = [0]*(len(fl))
		Data = []
		for ev = range(len(FitEvents)):
			DataEV = []
		  for ob = range(len(FitEvents[ev][1])):
				DataEv.append(fl.o[FitEvents[ev][1][ob]])
				#?????
			Data.append(numpy.matrix(DataEv))
		return Data
		
		def likelihood(self):
    return fitEKF.ekf(self.Data, self.M)
