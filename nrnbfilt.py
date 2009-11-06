from neuron import h
import noise
import models
import numpy
import fitEKF

class NrnBFilt(object):
  def __init__(self, ho):
    self.rf = ho
    ol = []
    P = noise.NoiseParams()
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
    h.cvode.states(s)
    assert(len(s) > 0)
    assert(len(vl) > 0)
    P.B = numpy.matrix(numpy.zeros((len(s), len(vl))))
    P.B[0,0] = .1
    P.InitialCov = numpy.eye(len(s))
    Obs = models.ObservationModel(P, 1000, ol)
    Sys = models.NeuronModel(P, 0, len(vl))
    Sys.Injection.erange(0.0, tlast, 1.0)
    self.M = models.Model(Sys, Obs, P)
    self.Data = self.__data(fl,self.M.FitEvents)
    print 'leave NrnBFilt'

  def __data(self,fl,FitEvents):
    counter = [1]*(len(fl))
    Data = []
    for elm in FitEvents:
      time = elm[0][-1]
      obindices = elm[1]
      DataEV = []
      for i in obindices:
        x = fl.o(i).xdat_
        y = fl.o(i).ydat_
        assert(time == x[counter[i]])
        DataEV.append(y[counter[i]])
        counter[i] += 1
      Data.append(numpy.matrix(DataEV))
    for i in range(len(fl)):
	assert(counter[i] == len(fl.o(i).xdat_))
    return Data
		
  def likelihood(self):
    print 'likelihood'
    x = fitEKF.ekf(self.Data, self.M)
    x = float(x[0,0])
    return x
