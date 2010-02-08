from neuron import h
import noise
import models
import numpy
import EKF
import math

class WrappedVal:
  def __init__(self, val):
    self.x = val

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
    P.B = numpy.matrix(numpy.zeros((len(s), len(s))))
    P.B[0,0] = 1
    self.processNoise = []
    for i in range(len(s)):
      self.processNoise.append(WrappedVal(P.B[i, i]))
    P.InitialCov = numpy.eye(len(s))
    Obs = models.ObservationModel(P, 1000, ol)
    Sys = models.NeuronModel(P, 0, len(vl))
    Sys.Injection.erange(0.0, tlast, 1.0)
    self.M = models.Model(Sys, Obs, P)
    self.Data = self.__data(fl,self.M)
    self.pf = self.getParmFitness()

  def __data(self,fl,M):
    counter = [0]*(len(fl))
    Data = []
    for idx, time in enumerate(M.collectionTimes):
      obindices = M.ObsNum[idx]
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
    print 'Collection Times\n', M.collectionTimes, '\nData\n', Data
    return Data

  def likelihood(self):
    x = EKF.ekf(self.Data, self.M)
    x = float(x)
    return -x

  def Etime(self):
    return h.Vector(EKF.Etime)

  def Ecenter(self, i):
    return h.Vector(EKF.Ecenter[int(i)])

  def Ewidth(self, i):
    return h.Vector(EKF.Ewidth[int(i)])

  def getParmFitness(self):
    # the ParmFitness instance that owns me.
    # there are probably not many so we can work forward from ParmFitness
    pfl = h.List('ParmFitness')
    for pf in pfl:
      for gi in pf.generatorlist:
        if gi.gen.hocobjptr() == self.rf.hocobjptr():
          return pf

  def getParm(self):
    #return Hoc Vector of current objective function parameters
    v = h.Vector()
    self.pf.doarg_get(v)
    return v

  def setParm(self, hvec):
    #assign current objective funtion parameters
    self.pf.parm(hvec)

  def fillPB(self, i):
    self.M.P.B[i,i] = self.processNoise[i].x
    print i, self.M.P.B

  def paramPanel(self):
    self.box = h.VBox()
    self.box.intercept(1)
    h.xpanel('')
    h.xlabel('Likelihood numerical parameters')
    h.xlabel('    Measurement noise')
    c =  self.M.Obs.C
    for o in c:
      h.xvalue('sigma: '+o.hpt.s(), (o, 'sigma'), 1)
    h.xlabel('    Process noise')
    s = h.Vector()
    h.cvode.states(s)
    sref = h.ref('')
    for i in range(len(s)):
      h.cvode.statename(i, sref, 1)
      h.xvalue('P.B[%d,%d]: '%(i,i) + sref[0], (self.processNoise[i], 'x'), 1, (self.fillPB, i))
    h.xpanel()
    self.box.intercept(0)
    self.box.map('Likelihood parameters')

