from neuron import h
import noise
import numpy
import EKF
import math
import sto
import detsys
import obs
import eve
import pickle
import cvodewrap
import fitglobals

class xstruct(object):
    def __init__(self):
        self.x = False

class WrappedVal:
    def __init__(self, val):
        self.x = val

class NrnBFilt(object):
    def dumpData(self):
        f = open('dataTrace.txt','w')
        Tlist = list(h.YFitness[0].xdat)
        Ylist = list(h.YFitness[0].ydat)
        for i in range(len(Tlist)):
            f.write('{0:f} {1:f}\n'.format(Tlist[i],Ylist[i]))
        f.close()

    def __init__(self, ho):
        h.mulfit_after_quad_pycallback = self.after_quad
        pc = h.ParallelContext()
        nhost = int(pc.nhost_bbs())
        if nhost > 1:
          fitglobals.verbose = 0
        self.xlvec = h.Vector()
        self.ylvec = h.Vector()
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
        cvodewrap.active(1)
        cvodewrap.states(s)
        assert(len(s) > 0)
        assert(len(vl) > 0)
        self.covGrowthTime = 100
        self.varTerm = 1
        self.processNoise = []
        self.Sdiag = []
        Sto = sto.StochasticModel(len(s),tlast)
        for i in range(len(s)):
            self.Sdiag.append(WrappedVal(Sto.InitialCovSqrt[i,i]))
        for i in range(len(s)):
            self.processNoise.append(WrappedVal(Sto.B[i, i]))
        Obs = obs.ObservationModel(ol)
        self.Eve = eve.EventTable(Sto,Obs)
        self.hhB = self.Eve.Sto.hhB
        self.nNa = self.Eve.Sto.nNa
        self.nK = self.Eve.Sto.nK
        self.Sys = detsys.NeuronModel()
        self.inj_invl = 1.0
        self.Eve.newInjectionInterval(self.inj_invl)
        # self.inj_invl_changed(Sys, P.tstop)
        # self.M = models.Model(Sys, Obs, P)
        self.Data = self.__data(fl,self.Eve)
        self.pf = self.getParmFitness()
        self.pf.verbose = fitglobals.verbose
        self.dlikedt = h.Vector()
        self.likefailed = False
        #CONSTRAINTS GUI INIT
        s = h.Vector()
        cvodewrap.states(s)
        nstates = len(s)
        self.geq0 = []
        self.leq1 = []
        self.sumto1 = []
        self.cvxopt_sel = True
        self.custom_sel = True
        self.nsums = 2
        for j in range(self.nsums):
            self.sumto1.append([])
        for i in range(nstates):
            self.geq0.append(xstruct())
            self.leq1.append(xstruct())
            for j in range(self.nsums):
                self.sumto1[j].append(xstruct())
        EKF.constraintsOn(self.geq0,self.leq1,self.sumto1)

    def __data(self,fl,Eve):
      n_trajectories = len(fl.o(0).ydat_)
      Data = []
      for i in range(n_trajectories):
        Data.append([])
      for i_trajec in range(n_trajectories):
        counter = [0]*(len(fl))
        for idx, time in enumerate(Eve.collectionTimes):
            obindices = Eve.ObsNum[idx]
            DataEV = []
            for i in obindices:
                x = fl.o(i).xdat_
                y = fl.o(i).ydat_[i_trajec]
                #print i, counter[i], time, x[counter[i]], time - x[counter[i]]
                assert(math.fabs(time - x[counter[i]]) < 1e-10)
                DataEV.append(y[counter[i]])
                counter[i] += 1
            Data[i_trajec].append(numpy.matrix(DataEV).T)
        for i in range(len(fl)):
            assert(counter[i] == len(fl.o(i).xdat_))
      #if fitglobals.verbose: print 'Collection Times\n', Eve.collectionTimes, '\nData\n', Data
      return Data

    def overwrite(self,Data):
        self.Data = Data

    def likelihood(self, trap_errors=False):
        self.ifchdat()
        # x = EKF.ekf(self.Data, self.Eve, self.Sys, DLikeDt_hvec = self.dlikedt)
        # x = float(x)
        # return -x
        x = 0
        if not trap_errors:
            for data in self.Data:
                x1 = EKF.ekf(data, self.Eve, self.Sys, DLikeDt_hvec = self.dlikedt)
                print x1
                x = x1
            x = float(x)
            #self.xlvec.append(self.getParm().x[0])
            #self.ylvec.append(x)
            return -x
        else:
            try:
                for data in self.Data:
                    x1 = EKF.ekf(data, self.Eve, self.Sys, DLikeDt_hvec = self.dlikedt)
                    print x1
                    x += x1
                x = float(x)
                #self.xlvec.append(self.getParm().x[0])
                #self.ylvec.append(x)
                return -x
            except:
                self.likefailed = True
                return float("1e12")

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

    def SUpper(self,i):
        v = h.Vector(EKF.Scenter[int(i)])
        return v.add(h.Vector(EKF.Swidth[int(i)]))

    def SLower(self,i):
        v = h.Vector(EKF.Scenter[int(i)])
        return v.sub(h.Vector(EKF.Swidth[int(i)]))

    def getParmFitness(self):
        # the ParmFitness instance that owns me.
        # there are probably not many so we can work forward from ParmFitness
        pfl = h.List('ParmFitness')
        for pf in pfl:
            for gi in pf.generatorlist:
                if gi.gen.hocobjptr() == self.rf.hocobjptr():
                    return pf

    def getParm(self):
        #return current objective function used parameters (log space)
        #return Hoc Vector
        v = h.Vector()
        self.pf.doarg_get(v)
        # print 'getParm', v[0]
        return v

    def getParmVal(self):
        #return current objective function used parameters (linear space)
        #return Hoc Vector
        # differs from getParm in that it is the values, not the log values.
        # differs from pf.argget, in that only the 'use' values are returned
        v = h.Vector()
        for o in self.pf.parmlist:
            if o.doarg > 0:
                v.append(o.val)
        return v

    def setParm(self, hvec):
        #assign current objective function used parameters (log space)
        # print 'setParm', hvec[0]
        self.pf.parm(hvec)

    def setParmVal(self, hvec):
        #assign current objective function used parameters (linear space)
        i = 0
        for o in self.pf.parmlist:
            if o.doarg > 0:
                o.val = hvec.x[i] 
                i += 1
                o.play_one()

    def fillS(self, i):
        print 'Sdiag', self.Sdiag[i].x
        self.Eve.Sto.InitialCovSqrt[i,i] = self.Sdiag[i].x
        print 'Cov Sqrt', self.Eve.Sto.InitialCovSqrt
        self.Initial_changed()
        print 'initial cov', i, self.Eve.Sto.InitialCov

    def fillPB(self, i):
        if i == -1:
            for j in range(self.Eve.Sto.B.shape[0]):
                self.Eve.Sto.B[j,j] = self.processNoise[j].x
        else:
            self.Eve.Sto.B[i,i] = self.processNoise[i].x
        print i, self.Eve.Sto.B

    def inj_invl_changed(self):
        self.Eve.newInjectionInterval(self.inj_invl)
        self.Data = self.__data(self.rf.fitnesslist,self.Eve)

    def Initial_changed(self):
        self.Eve.Sto.updateInitial()

    def constraintsPanel(self):
        self.box = h.HBox()
        self.box.intercept(1)
        self.box.ref(self)
        s = h.Vector()
        cvodewrap.states(s)
        nstates = len(s)
        sref = h.ref('')
        h.xpanel("")
        h.xlabel('0<=')
        for i in range(nstates):
            h.xcheckbox('',(self.geq0[i],'x'), self.constraintsButton)
        h.xpanel()
        h.xpanel("")
        h.xlabel('>=1')
        for i in range(nstates):
            cvodewrap.statename(i,sref,1)
            h.xcheckbox(sref[0],(self.leq1[i],'x'),self.constraintsButton)
        h.xpanel()
        for j in range(self.nsums):
            h.xpanel("")
            h.xlabel("S%d"%j)
            for i in range(nstates):
                h.xcheckbox('',(self.sumto1[j][i],'x'),self.constraintsButton)
            h.xpanel()
        h.xpanel("")
        h.xbutton("Add Sum-to-1", self.constraintsButton)
        h.xbutton("Remove Empty", self.constraintsButton)
        h.xbutton("Close", self.constraintsButton)
        h.xlabel('QP Solver:')
        h.xstatebutton('cvxopt',(self,'cvxopt_sel'),self.constraintsButton)
        h.xstatebutton('custom',(self,'custom_sel'),self.constraintsButton)
        h.xpanel()
        self.box.intercept(0)
        self.box.map("Constraints")

    def constraintsButton(self):
            EKF.constraintsOn(self.geq0,self.leq1,self.sumto1)

    def hhBButton(self):
        self.Eve.Sto.hhB = self.hhB
        self.Eve.Sto.nNa = self.nNa
        self.Eve.Sto.nK = self.nK
        if self.Eve.Sto.hhB:
            print 'Fox & Lu ON: nNa', self.Eve.Sto.nNa, 'nK', self.Eve.Sto.nK
        else:
            print 'Fox & Lu OFF'
        self.fillPB(-1)

    def inc_nsums(self):
        self.nsums += 1
        s = h.Vector()
        cvodewrap.states(s)
        nstates = len(s)
        new = []
        for i in range(nstates):
            new.append(xstruct())
        self.sumto1.append(new)
        print 'nsums', self.nsums
        print 'sumto1',self.sumto1


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
        cvodewrap.states(s)
        sref = h.ref('')
        for i in range(len(s)):
            cvodewrap.statename(i, sref, 1)
            h.xvalue('Diffusion Coeff[%d,%d]: '%(i,i) + sref[0], (self.processNoise[i], 'x'), 1, (self.fillPB, i))
        h.xcheckbox('Fox & Lu Diffusion (for Hodgkin-Huxley)?',(self,'hhB'), self.hhBButton)
        h.xvalue('  Fox & Lu: Number Na Channels', (self,'nNa'), 1, self.hhBButton)
        h.xvalue('  Fox & Lu: Number K Channels', (self,'nK'), 1, self.hhBButton)
        h.xlabel('    Initial Uncertainty')
        for i in range(len(s)):
            print i
            cvodewrap.statename(i, sref, 1)
            h.xvalue('Initial Stand Dev[%d]: '%i + sref[0], (self.Sdiag[i],'x'), 1, (self.fillS,i))
        h.xbutton('Show state funnels', self.show_state_funnels)
        h.xpanel()
        self.box.intercept(0)
        self.box.map('Likelihood parameters')

    def data_change(self):
        print "data_change", self
        inj_invl = self.inj_invl
        covGrowthTime = self.covGrowthTime
        varTerm = self.varTerm
        hhB = self.hhB
        nNa = self.nNa
        nK = self.nK
        c = self.Eve.Obs.C
        scl = self.Eve.Sto.scale
        pn = self.processNoise
        sd = self.Sdiag
        g = self.g

        self.__init__(self.rf)

        self.inj_invl = inj_invl
        self.covGrowthTime = covGrowthTime
        self.varTerm = varTerm
        self.hhB = hhB
        self.nNa = nNa
        self.nK = nK
        self.Eve.Sto.hhB = hhB
        self.Eve.Sto.nNa = nNa
        self.Eve.Sto.nK = nK
        self.Eve.Sto.scale = scl
        cnew = self.Eve.Obs.C
        for i in range(len(c)):
            cnew[i].sigma = c[i].sigma
        for i in range(len(pn)):
            self.processNoise[i].x = pn[i].x
            self.Eve.Sto.B[i,i] = pn[i].x
        for i in range(len(sd)):
            self.Sdiag[i].x = sd[i].x
            self.Eve.Sto.InitialCovSqrt[i,i] = sd[i].x
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
            #self.Scenter(i).line(g[i], t)
            self.SUpper(i).line(g[i], t)
            self.SLower(i).line(g[i], t)
            if self.g == None:
                g[i].exec_menu('View = plot')
        self.g = g

    def save_session(self, nameprefix):
        f = open(nameprefix + '.lkl', 'w')
        # f.write('here is some info from %s\n' % self.save_session)
        pickle.dump(self.geq0,f)
        pickle.dump(self.leq1,f)
        pickle.dump(self.sumto1,f)
        pickle.dump(self.cvxopt_sel,f)
        pickle.dump(self.custom_sel,f)
        pickle.dump(self.nsums,f)
        c = self.Eve.Obs.C
        for o in c:
            pickle.dump(o.sigma,f)
        pickle.dump(self.inj_invl,f)
        pickle.dump(self.processNoise,f)
        pickle.dump(self.covGrowthTime,f)
        pickle.dump(self.varTerm,f)
        pickle.dump(self.Sdiag,f)
        pickle.dump(self.Eve.Sto.hhB,f)

    def restore_session(self, nameprefix):
        try:
            f = open(nameprefix + '.lkl', 'r')
        except:
            return
        # for line in f:
        #    print 'restored: ', line
        self.geq0 = pickle.load(f)
        self.leq1 = pickle.load(f)
        self.sumto1 = pickle.load(f)
        self.cvxopt_sel = pickle.load(f)
        self.custom_sel = pickle.load(f)
        self.nsums = pickle.load(f)
        c = self.Eve.Obs.C
        for o in c:
            o.sigma = pickle.load(f)
        self.inj_invl = pickle.load(f)
        self.processNoise = pickle.load(f)
        self.covGrowthTime = pickle.load(f)
        self.varTerm = pickle.load(f)
        try:
            self.Sdiag = pickle.load(f)
        except:
            print 'Initial Standard Deviations not yet saved'
        try:
            self.Eve.Sto.hhB = pickle.load(f)
        except:
            print 'HH Button Flag not yet saved'

    def after_quad(self, n):
        #print 'after_quad ', n
        pass
