from neuron import h
import noise
import numpy
import EKF
import math
import sto
import sys
import obs
import eve

class WrappedVal:
  def __init__(self, val):
    self.x = val

class NrnBFilt(object):
  def __init__(self, ho):
    self.rf = ho
    ol = []
    vl = self.rf.yvarlist
    fl = self.rf.fitnesslist
    tlast = 0
    for i in range(len(vl)):
      tl = list(fl.o(i).xdat_)
      o = obs.NeuronObservable(vl.o(i), tl)
      o.sigma = 0.01
      ol.append(o)
      if (tlast < tl[-1]):
        tlast = tl[-1]
    s = h.Vector()
    h.cvode.active(1)
    h.cvode.states(s)
    assert(len(s) > 0)
    assert(len(vl) > 0)
    Sto = sto.StochasticModel(len(s),tlast)
    self.processNoise = []
    for i in range(len(s)):
      self.processNoise.append(WrappedVal(Sto.B[i, i]))
    Obs = obs.ObservationModel(ol)
    self.Eve = eve.EventTable(Sto,Obs)
    self.Sys = sys.NeuronModel()
    # self.inj_invl = 1.0
    # self.inj_invl_changed(Sys, P.tstop)
    # self.M = models.Model(Sys, Obs, P)
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
    x = EKF.ekf(self.Data, self.Eve, self.Sys)
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

  def inj_invl_changed(self, sys, tstop):
    sys.Injection.erange(0.0, tstop, self.inj_invl)

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
    h.xvalue('Injection interval', (self, 'inj_invl'), 1, (self.Eve.newInjectionInterval, self.inj_invl))
    s = h.Vector()
    h.cvode.states(s)
    sref = h.ref('')
    for i in range(len(s)):
      h.cvode.statename(i, sref, 1)
      h.xvalue('P.B[%d,%d]: '%(i,i) + sref[0], (self.processNoise[i], 'x'), 1, (self.fillPB, i))
    h.xpanel()
    self.box.intercept(0)
    self.box.map('Likelihood parameters')

