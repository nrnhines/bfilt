import pylab
import numpy
import random
import math
import scipy
import scipy.linalg
import numpy.linalg

def ch3Q(alpha01, beta01, alpha12, beta12):
    Q = numpy.matrix([[-alpha01, alpha01, 0], [beta01, -(beta01 + alpha12), alpha12], [0, beta12, -beta12]])
    return Q

def ch3Qv(V, tau01=2.,tau12=4.,Vhalf01=-20.,Vhalf12=-25.,Vslope01=1.,Vslope12=1.):
    inf01 = 1./(1. + math.exp((1./Vslope01)*(Vhalf01 - V)))
    inf12 = 1./(1. + math.exp((1./Vslope12)*(Vhalf12 - V)))
    alpha01 = inf01/tau01
    beta01 = (1./tau01)-alpha01
    alpha12 = inf12/tau12
    beta12 = (1./tau12)-alpha12
    Q = ch3Q(alpha01, beta01, alpha12, beta12)
    # print Q
    return Q

def ch3hmm(V0=-65.,V1=20.,tau01=2.,tau12=4.,Vhalf01=-20.,Vhalf12=-25.,Vslope01=1.,Vslope12=1.,sigma=0.001):
    if V0 < -40:
        pstates = [1.0, 0.0, 0.0]  # Saves an expensive eig computation
    elif V0 > 0:
        pstates = [0.0, 0.0, 1.0]  # Saves an expensive eig computation
    else:
        Q0 = ch3Qv(V0, tau01, tau12,Vhalf01,Vhalf12,Vslope01,Vslope12)
        pstates = equilibrium(Q0)
    output = [0.0, 0.0, 1.0]
    Q = ch3Qv(V1, tau01, tau12)
    H = HMM(pstates,output,Q,sigma)
    return H

def equilibrium(Q):
    (V,D) = numpy.linalg.eig(Q.T)
    # print 'Q.T',Q.T
    # print 'V', V
    # print 'D', D
    # print 'd', D[:,0]
    # print 'Q.T*d', Q.T*D[:,0]
    # Find index of eigenvalue 0
    m = 1.0
    for i in range(V.shape[0]):
        if math.fabs(V[i]) < m:
            mi = i
            m = math.fabs(V[i])
    eigvect = D[:,mi]
    # Find index of eigenvector component with largest magnitude (to determine sign)
    m = 1e-7
    for i in range(V.shape[0]):
        if math.fabs(eigvect[i,0]) > m:
            mi = i
            m = math.fabs(eigvect[i,0])
    # Reverse sign if necessary.  (if so all components should be negative)
    if eigvect[mi,0] < 0:
        eigvect = -eigvect
    pstates = []
    for i in range(V.shape[0]):
        pstates.append(eigvect[i,0])
    normalization = sum(pstates)
    for i in range(V.shape[0]):
        pstates[i] = pstates[i]/normalization
    tol1 = 1e-12
    tol = 1e-6
    # print 'pstates', pstates
    # print 'Q0', Q
    # print 'pstates*exp(Q)', pstates*scipy.linalg.expm(Q)
    for i in range(len(pstates)):
        assert(pstates[i].imag == 0)
        # assert(math.fabs(pstates[i].imag)<tol1)
        # pstates[i] = pstates[i].real
        # print "HERE: pstates[i]", pstates[i]
        assert(pstates[i] >= -tol1)
        assert(pstates[i] <= 1.0+tol1)
        if pstates[i] < 0:
            pstates[i] = 0
        if pstates[i] > 1:
            pstates[i] = 1
    assert(sum(pstates)<1.0+tol)
    assert(sum(pstates)>1.0-tol)
    return pstates

class HMM(object):
    def __init__(self, pstates, output, Q, sigma=0.001):
        self.simmed = False
        self.liked = False
        self.nstates = len(pstates)
        self.Q = Q
        self.pstates = pstates
        self.init = numpy.matrix(pstates)
        self.output = output
        self.dt = None
        self.skip = None
        self.sigma = sigma #std of measurement error
        self.simStates = None
        self.simOut = None
        self.nsamples = None
        self.simData = None
        self.fitData = None
        self.center = None
        self.width = None
        self.R = random.Random()
        tol = 1e-7
        assert self.init.shape[1] == self.Q.shape[0]
        assert self.init.shape[1] == self.Q.shape[1]
        assert self.init.shape[0] == 1
        assert len(output) == self.init.shape[1]
        assert math.fabs(sum(pstates) - 1.0) < tol

    def select(self,mat,row):
        p = self.R.random()
        rowsum = 0
        for i in range(mat.shape[1]):
            rowsum += mat[row,i]
            if p < rowsum:
                return i
        assert False

    def sim(self, seeds=[0], dt=0.1, tstop=20):
        if not type(seeds) is list:
            seeds = [seeds]
        self.simseeds = seeds
        self.ntraj = len(self.simseeds)
        self.simdt = dt
        self.simtstop = tstop
        self.simtrans = scipy.linalg.expm(dt*self.Q)
        tol = 1e-7
        for i in range(self.nstates):
            rowsum = 0
            for j in range(self.nstates):
                rowsum += self.simtrans[i,j]
            assert math.fabs(rowsum - 1.0) < tol
        assert(tstop>0.0)
        print "Checks out..."
        self.nsamples = int(math.ceil(tstop/dt))
        self.simData = []
        for seed in self.simseeds:
            self.R.seed(seed)
            self.simStates = []
            self.simStates.append(self.select(self.init,0))
            for i in range(self.nsamples-1):
                self.simStates.append(self.select(self.simtrans,self.simStates[-1]))
            # print 'States as a function of discrete time:', self.simStates
            self.simOut = []
            for s in self.simStates:
                self.simOut.append(self.output[s])
            # print 'Output with out noise:', self.simOut
            simDataX = []
            simDataT = []
            t = 0.0
            for o in self.simOut:
                t+=dt
                simDataT.append(t)
                simDataX.append(o + self.R.normalvariate(0,self.sigma))
            # print 'simData:', self.simData
            self.simData.append((simDataT,simDataX))
        self.simmed = True

    def normpdf(self,x,m,sigma):
        return (1/(math.sqrt(2*math.pi)*sigma))*math.exp(-(x-m)*(x-m)/(2*sigma*sigma))

    def initializeErrorBars(self):
        self.center = []
        self.width = []
        self.time = [0.0]

    def saveErrorBars(self, pmf, t):
        E = 0.0
        for i in range(len(self.output)):
            E += self.output[i]*pmf[0,i]
        V = 0.0
        for i in range(len(self.output)):
            V += ((self.output[i] - E)**2)*pmf[0,i]
        self.center.append(E)
        self.width.append(math.sqrt(V+self.sigma**2))
        if t>0.0:  # self.time initialized as [0.0]
            self.time.append(t)

    def predict(self,skips,extraskipdt,extratrans,inter):
        # print 'time', self.time[-1], 'pmf', inter
        self.saveErrorBars(inter,self.time[-1])
        if not self.plotskipdt == None:
            for i in range(skip):
                inter = inter*self.skiptrans
                self.saveErrorBars(inter,self.time[-1]+self.plotskipdt)
        inter = inter*extratrans
        self.saveErrorBars(inter,self.time[-1]+extraskipdt)
        return inter

    def update(self,datapoint,prior):
        weight = []
        for o in self.output:
            weight.append(self.normpdf(datapoint - o,0,self.sigma))
        new = []
        for i in range(len(weight)):
            new.append(weight[i]*prior[0,i])
        marg = sum(new)
        posterior = numpy.matrix(new)/marg
        # print 'datapoint', datapoint
        # print 'weight', weight
        # print 'prior', prior
        # print 'posterior', posterior
        return (posterior, marg)

    def likelihood(self, fitData, plotskipdt=None):
        self.fitData = fitData
        total = 0
        for fD in fitData:
            ts = fD[0]
            xs = fD[1]
            assert(len(ts) == len(xs))
            nData = len(xs)
            self.processTimeIntervals(ts, plotskipdt)
            self.initializeErrorBars()
            pmf = self.init
            sll = 0
            for i in range(nData):
                # print i, i*self.dt #, self.fitData[i]
                pre = self.predict(self.plotskips[i],self.extraskipdt[i],self.extratrans[i],pmf)
                (pmf, lk) = self.update(xs[i], pre)
                sll += math.log(lk)
            total += sll
        self.liked = True
        return total

    def processTimeIntervals(self, ts, plotskipdt):
        # Still need to make first time = 0 possible
        self.plotskipdt = plotskipdt
        nData = len(ts)
        tpre = 0.0
        self.dts = []
        for t in ts:
            newdt = t-tpre
            assert(newdt>0) #Must make first time OK for t0 = 0 newdt=0
            self.dts.append(newdt)
            tpre = t
        mindt = min(self.dts)
        maxdt = max(self.dts)
        dttol = 1e-7
        skiptol = 1e-7
        if (not self.plotskipdt == None) and self.plotskipdt < maxdt - dttol:  # I.e. we want to plot likelihood at higher sample freq
            self.plotskiptrans = scipy.linalg.expm(plotskipdt*self.Q)
            self.plotskips = []
            self.extraskipdt = []
            self.extratrans = []
            for dt in dts:
                nskip = int(math.floor(dt/plotskipdt))
                assert(nskip > -1)
                self.extradt.append(dt - nskip*plotskipdt)
                assert(self.extraskipdt[-1] > - skiptol)
                self.plotskips.append(nskip)
                if self.extraskipdt[-1] > skiptol:
                    self.extratrans.append(scipy.linalg.expm(self.extraskipdt[-1]*self.Q))
                else:
                    self.extratrans = None
        elif max(self.dts) - min(self.dts) < dttol:  # i.e. if all dts are approximately the same
            self.plotskipdt = None
            self.extratrans = [scipy.linalg.expm(self.dts[-1]*self.Q)]*nData
            self.extraskipdt = [self.dts[-1]]*nData
            self.plotskips = [0]*nData
        else:
            self.plotskipdt = None
            self.plotskips = []
            self.extraskipdt = []
            self.extratrans = []
            for dt in dts:
                self.extratrans.append(scipy.linalg.expm(dt*self.Q))
                self.extraskipdt.append(dt)
                self.plotskips.append(0)

    def simplot(self,num=0):
        assert(self.simmed)
        x = numpy.array(self.simData[num][0])
        y = numpy.array(self.simData[num][1])
        pylab.hold(False)
        pylab.plot(x,y)

    def plot(self):
        #MIGHT NOT WORK AFTER MULTIPLE TRAJS ADDED
        assert(self.liked)
        x = numpy.array(self.fitData[0])
        y = numpy.array(self.fitData[1])
        pylab.hold(False)
        pylab.plot(x,y)
        xf = self.time
        yfh = []
        yfl = []
        for i in range(len(self.center)):
            yfh.append(self.center[i] + self.width[i])
            yfl.append(self.center[i] - self.width[i])
        pylab.hold(True)
        pylab.plot(xf,yfh,'g')
        pylab.plot(xf,yfl,'g')
