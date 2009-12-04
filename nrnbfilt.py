from neuron import h
import noise
import models
import numpy
import fitEKF
import math

class NrnBFilt(object):
  def __init__(self, ho):
    self.rf = ho
    ol = []
    P = noise.NoiseParams()
    P.dt = 0.0
    vl = self.rf.yvarlist
    fl = self.rf.fitnesslist
    tlast = 0
    for i in range(len(vl)):
      tl = list(fl.o(i).xdat_)
      o = models.NeuronObservable(vl.o(i), P, i, tl)
      o.sigma = 0.01
      ol.append(o)
      if (tlast < tl[-1]):
        tlast = tl[-1]
    P.tstop = tlast
    s = h.Vector()
    h.cvode.active(1)
    h.cvode.states(s)
    assert(len(s) > 0)
    assert(len(vl) > 0)
    P.B = numpy.matrix(numpy.zeros((len(s), len(vl))))
    P.B[0,0] = .01
    P.InitialCov = numpy.eye(len(s))
    Obs = models.ObservationModel(P, 1000, ol)
    Sys = models.NeuronModel(P, 0, len(vl))
    Sys.Injection.erange(0.0, tlast, 1.0)
    self.M = models.Model(Sys, Obs, P)
    self.Data = self.__data(fl,self.M.FitEvents)
    #print 'leave NrnBFilt'

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
        #print i, counter[i], time, x[counter[i]], time - x[counter[i]]
        assert(math.fabs(time - x[counter[i]]) < 1e-10)
        DataEV.append(y[counter[i]])
        counter[i] += 1
      Data.append(numpy.matrix(DataEV).T)
    for i in range(len(fl)):
	assert(counter[i] == len(fl.o(i).xdat_))
    print 'FitEvents\n', FitEvents, '\nData\n', Data
    return Data
		
  def likelihood(self):
    x = fitEKF.ekf(self.Data, self.M)
    x = float(x)
    return -x
