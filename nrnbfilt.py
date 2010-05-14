from neuron import h
import noise
import numpy
import EKF
import math
import sto
import detsys
import obs
import eve

class WrappedVal:
    def __init__(self, val):
        self.x = val

class NrnBFilt(object):
    def __init__(self, ho):
        self.g = None
        self.rf = ho
        ol = []
        vl = self.rf.yvarlist
        fl = self.rf.fitnesslist
        tlast = 0
        self.n_chdat = 0
        for i in range(len(vl)):
            self.n_chdat += fl.o(i).n_chdat
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
        self.covGrowthTime = 100
        self.varTerm = 1
        Sto = sto.StochasticModel(len(s),tlast,self.covGrowthTime,self.varTerm)
        self.processNoise = []
        for i in range(len(s)):
            self.processNoise.append(WrappedVal(Sto.B[i, i]))
        Obs = obs.ObservationModel(ol)
        self.Eve = eve.EventTable(Sto,Obs)
        self.Sys = detsys.NeuronModel()
        self.inj_invl = 1.0
        self.Eve.newInjectionInterval(self.inj_invl)
        # self.inj_invl_changed(Sys, P.tstop)
        # self.M = models.Model(Sys, Obs, P)
        self.Data = self.__data(fl,self.Eve)
        self.pf = self.getParmFitness()
        self.dlikedt = h.Vector()
    
    def __data(self,fl,Eve):
        counter = [0]*(len(fl))
        Data = []
        for idx, time in enumerate(Eve.collectionTimes):
            obindices = Eve.ObsNum[idx]
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
        print 'Collection Times\n', Eve.collectionTimes, '\nData\n', Data
        return Data
    
    def likelihood(self):
        self.ifchdat()
        x = EKF.ekf(self.Data, self.Eve, self.Sys, DLikeDt_hvec = self.dlikedt)
        x = float(x)
        return -x
    
    def ifchdat(self):
        fl = self.rf.fitnesslist
        n = 0
        for i in range(len(fl)):
            n += fl.o(i).n_chdat
        if n != self.n_chdat:
            self.data_change()

    def Etime(self):
        return h.Vector(EKF.Etime)
    
    def Ecenter(self, i):
        return h.Vector(EKF.Ecenter[int(i)])
    
    def Ewidth(self, i):
        return h.Vector(EKF.Ewidth[int(i)])
    
    def Scenter(self, i):
        return h.Vector(EKF.Scenter[int(i)])
    
    def Swidth(self, i):
        return h.Vector(EKF.Swidth[int(i)])
    
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
        self.Eve.Sto.B[i,i] = self.processNoise[i].x
        self.Initial_changed()
        print i, self.Eve.Sto.B
    
    def inj_invl_changed(self):
        self.Eve.newInjectionInterval(self.inj_invl)
        self.Data = self.__data(self.rf.fitnesslist,self.Eve)
    
    def Initial_changed(self):
        self.Eve.Sto.updateInitial(self.covGrowthTime,self.varTerm)
    
    def paramPanel(self):
        self.box = h.VBox()
        self.box.intercept(1)
        h.xpanel('')
        h.xlabel('Likelihood numerical parameters')
        h.xlabel('    Measurement noise')
        c =  self.Eve.Obs.C
        for o in c:
            h.xvalue('sigma: '+o.hpt.s(), (o, 'sigma'), 1)
        h.xlabel('    Process noise')
        h.xvalue('Injection interval', (self, 'inj_invl'), 1, self.inj_invl_changed)
        s = h.Vector()
        h.cvode.states(s)
        sref = h.ref('')
        for i in range(len(s)):
            h.cvode.statename(i, sref, 1)
            h.xvalue('P.B[%d,%d]: '%(i,i) + sref[0], (self.processNoise[i], 'x'), 1, (self.fillPB, i))
        h.xlabel('    Initial Uncertainty')
        h.xvalue('Covariance Growth Time', (self, 'covGrowthTime'), 1, self.Initial_changed)
        h.xvalue('Additional Variance', (self, 'varTerm'), 1, self.Initial_changed)
        h.xbutton('Show state funnels', self.show_state_funnels)
        h.xpanel()
        self.box.intercept(0)
        self.box.map('Likelihood parameters')

    def data_change(self):
        inj_invl = self.inj_invl
        covGrowthTime = self.covGrowthTime
        varTerm = self.varTerm
        c = self.Eve.Obs.C
        pn = self.processNoise
        g = self.g

        self.__init__(self.rf)

        self.inj_invl = inj_invl
        self.covGrowthTime = covGrowthTime
        self.varTerm = varTerm
        cnew = self.Eve.Obs.C
        for i in range(len(c)):
            cnew[i].sigma = c[i].sigma
        for i in range(len(pn)):
            self.processNoise[i].x = pn[i].x
            self.Eve.Sto.B[i,i] = pn[i].x
        self.Initial_changed()
        self.inj_invl_changed()       
        self.g = g

    def show_state_funnels(self):
        t = self.Etime()
        if self.g == None:
            g = []
            for i in range(len(EKF.Scenter)):
                g.append(h.Graph())
        else:
            g = self.g
        for i in range(len(EKF.Scenter)):
            self.Scenter(i).line(g[i], t)
            if self.g == None:
                g[i].exec_menu('View = plot')
        self.g = g

    def save_session(self, nameprefix):
	f = open(nameprefix + '.lkl', 'w')
        f.write('here is some info from %s\n' % self.save_session)

    def restore_session(self, nameprefix):
	f = open(nameprefix + '.lkl', 'r')
	for line in f:
            print 'restored: ', line
