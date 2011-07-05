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

def ch3Qv(V,tau01=2.,tau12=4.,Vhalf01=-20.,Vhalf12=-25.,Vchar01=1.,Vchar12=1.):
    # print "V", V
    inf01 = 1./(1. + math.exp((1./Vchar01)*(Vhalf01 - V)))
    inf12 = 1./(1. + math.exp((1./Vchar12)*(Vhalf12 - V)))
    alpha01 = inf01/tau01
    beta01 = (1./tau01)-alpha01
    alpha12 = inf12/tau12
    beta12 = (1./tau12)-alpha12
    Q = ch3Q(alpha01, beta01, alpha12, beta12)
    # print 'Q', Q
    return Q

def ch3hmm(V0=-65.,V1=20.,tau01=2.,tau12=4.,Vhalf01=-20.,Vhalf12=-25.,Vchar01=1.,Vchar12=1.,sigma=0.001):
    if V0 < min(Vhalf01,Vhalf12) - 15*max(abs(Vchar01),abs(Vchar12)):
        pstates = [1.0, 0.0, 0.0]  # Saves an expensive eig computation
    elif V0 > max(Vhalf01,Vhalf12) + 15*max(abs(Vchar01),abs(Vchar12)):
        pstates = [0.0, 0.0, 1.0]  # Saves an expensive eig computation
    else:
        Q0 = ch3Qv(V0,tau01,tau12,Vhalf01,Vhalf12,Vchar01,Vchar12)
        pstates = equilibrium(Q0)
    # print "initial state", pstates
    output = [0.0, 0.0, 1.0]
    Q = ch3Qv(V1,tau01,tau12,Vhalf01,Vhalf12,Vchar01,Vchar12)
    H = HMM(pstates,output,Q,sigma)
    return H

def ch3chain(V0,V1,tau01=2.,tau12=4.,Vhalf01=-20,Vhalf12=-25,Vchar01=1,Vchar12=1.,sigma=0.001):
    if V0 < min(Vhalf01,Vhalf12) - 15*max(abs(Vchar01),abs(Vchar12)):
        pstates = [1.0, 0.0, 0.0]  # Saves an expensive eig computation
    elif V0 > max(Vhalf01,Vhalf12) + 15*max(abs(Vchar01),abs(Vchar12)):
        pstates = [0.0, 0.0, 1.0]  # Saves an expensive eig computation
    else:
        Q0 = ch3Qv(V0,tau01,tau12,Vhalf01,Vhalf12,Vchar01,Vchar12)
        pstates = equilibrium(Q0)
    # print "Equilibrium Q0", pstates
    output = [0.0, 0.0, 1.0]
    Q = []
    for V in V1:
        Q.append(ch3Qv(V,tau01,tau12,Vhalf01,Vhalf12,Vchar01,Vchar12))
        H = HMMChain(pstates,output,Q,sigma)
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

class HMMChain(object):
    def __init__(self, pstates, output, Q, sigma=0.001):
        self.simmed = False
        self.liked = False
        self.pstates = pstates
        self.output = output
        self.Q = Q
        self.sigma = sigma
        self.nstates = len(pstates)
        self.R = random.Random()
        
    def sim(self, seeds=[0], dt=0.1, tstops=None):
        if tstops == None:
            tstops = [20.]*len(self.Q)
        assert len(tstops) > 0
        self.HMMLinks = []
        # Put links together
        for seedi in range(len(seeds)):
            self.HMMLinks.append([])
            for stopj in range(len(tstops)):
                # print "Q#", stopj%len(self.Q)
                self.HMMLinks[-1].append(HMM(self.pstates,self.output,self.Q[stopj%len(self.Q)],self.sigma,self.R))
        if len(tstops) < len(self.Q):
            print "Warning: fewer stop-times than Q's, truncating protocol"
        if len(tstops) > len(self.Q):
            print "Warning: more stop-times than Q's, looping protocol"
        if not type(seeds) is list:
            seeds = [seeds]
        self.simseeds = seeds
        self.ntraj = len(self.simseeds)
        self.simdt = dt
        self.simtstops = tstops
        for j in range(len(self.simseeds)):
            self.R.seed(self.simseeds[j])
            nextFirstState = None
            for i in range(len(tstops)):
                # QIndex = i%len(self.Q)  # QIndex = i mod len(Q)   
                self.HMMLinks[j][i].sim(None,self.simdt,[self.simtstops[i]],nextFirstState)
                nextFirstState = self.HMMLinks[j][i].simStates[-1]
        self.simmed = True
    
    def likelihood(self, fitChain, plotskipdt=None):
        self.fitChain = fitChain
        j = 0
        total = 0
        for j in range(len(fitChain.HMMLinks)):
            nextinitpmf=None
            for i in range(len(fitChain.HMMLinks[j])):
                total+=self.HMMLinks[j][i].likelihood(fitChain.HMMLinks[j][i].simData,plotskipdt,nextinitpmf)
                nextinitpmf = self.HMMLinks[j][i].likefinalpmf
        return total
    
    def simplot(self,num=0):
        assert(self.simmed)
        seednum = 0
        hold = False
        for oneseed in self.HMMLinks: 
            if seednum == num:
                redSeed = oneseed
            else:
                tOrigin = 0
                for H in oneseed:
                    col = 'b'
                    H.simplot(0,tOrigin,col,col,hold)
                    hold = True
                    tOrigin += H.simData[0][0][-1]
            seednum += 1
        tOrigin = 0.
        for H in redSeed:
            col = 'r'
            H.simplot(0,tOrigin,col,col,hold)
            hold = True
            tOrigin += H.simData[0][0][-1]

class HMM(object):
    def __init__(self, pstates, output, Q, sigma=0.001, R=None):
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
        if R == None:
            self.R = random.Random()
        else:
            self.R = R
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

    def sim(self, seeds=[0], dt=0.1, tstops=[20], firstState=None):
        if type(tstops) is list:
            tstop = tstops[0]
        else:
            tstop = tstops
        if not seeds == None:
            if not type(seeds) is list:
                seeds = [seeds]
            self.simseeds = seeds
            self.ntraj = len(self.simseeds)
        else:
            self.simseeds = [None]
            self.ntraj = 1
        self.simdt = dt
        self.firstState = firstState
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
            if not seed == None:
                self.R.seed(seed)
            self.simStates = []
            if self.firstState == None:
                # Pick a random state according to pmf self.init=numpy.matrix(pstates)
                self.simStates.append(self.select(self.init,0))
            else:
                self.simStates.append(self.firstState)
            for i in range(self.nsamples-1):
                # pick a random state as a function of pmf simState[-1]st row of self.simtrans
                self.simStates.append(self.select(self.simtrans,self.simStates[-1]))
            # print 'States as a function of discrete time:', self.simStates
            self.simOut = []
            for s in self.simStates:
                self.simOut.append(self.output[s])
            # print 'Output without noise:', self.simOut
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
        assert (inter>=0).all()
        self.saveErrorBars(inter,self.time[-1])
        if not self.plotskipdt == None:
            for i in range(skip):
                inter = inter*self.skiptrans
                assert (inter>=0).all()
                self.saveErrorBars(inter,self.time[-1]+self.plotskipdt)
        assert (inter*extratrans >= 0).all()
        inter = inter*extratrans
        assert (inter>=0).all()
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
        assert marg>0
        return (posterior, marg)

    def likelihood(self, fitData, plotskipdt=None, initpmf=None):
        self.fitData = fitData
        total = 0
        # For chains fitData is a list with one element
        for fD in fitData:
            ts = fD[0]
            xs = fD[1]
            assert(len(ts) == len(xs))
            nData = len(xs)
            self.processTimeIntervals(ts, plotskipdt)
            self.initializeErrorBars()
            if initpmf == None:
                pmf = self.init
            else:
                pmf = initpmf
            sll = 0
            for i in range(nData):
                # print i, i*self.dt #, self.fitData[i]
                pre = self.predict(self.plotskips[i],self.extraskipdt[i],self.extratrans[i],pmf)
                (pmf, lk) = self.update(xs[i], pre)
                assert (pmf>=0).all()
                sll += math.log(lk)
            total += sll
        self.likefinalpmf=pmf
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
        # INSIDE THIS IF STATEMENT NOT TESTED, ELIF YES, ELSE NO
        if (not self.plotskipdt == None) and self.plotskipdt < maxdt - dttol:  # I.e. we want to plot likelihood at higher sample freq
            print "Case 0"
            self.plotskiptrans = scipy.linalg.expm(plotskipdt*self.Q)
            self.plotskips = []
            self.extraskipdt = []
            self.extratrans = []
            for dt in dts:
                nskip = int(math.floor(dt/plotskipdt))
                assert(nskip > -1)
                self.extraskipdt.append(dt - nskip*plotskipdt) #extraSKIPdt  (SKIP added later)
                assert(self.extraskipdt[-1] > - skiptol)
                self.plotskips.append(nskip)
                if self.extraskipdt[-1] > skiptol:
                    self.extratrans.append(scipy.linalg.expm(self.extraskipdt[-1]*self.Q))
                else:
                    self.extratrans = None
        elif max(self.dts) - min(self.dts) < dttol:  # i.e. if all dts are approximately the same
            print "Case 1"
            self.plotskipdt = None
            self.extratrans = [scipy.linalg.expm(self.dts[-1]*self.Q)]*nData  #makes a sequence nData long
            self.extraskipdt = [self.dts[-1]]*nData
            self.plotskips = [0]*nData
        else:
            print "Case 2"
            self.plotskipdt = None
            self.plotskips = []
            self.extraskipdt = []
            self.extratrans = []
            for dt in dts:
                self.extratrans.append(scipy.linalg.expm(dt*self.Q))
                self.extraskipdt.append(dt)
                self.plotskips.append(0)

    def simplot(self,num=0,tOrigin=0.,colY='r',colN='b',hold=False):
        assert(self.simmed)
        pylab.hold(hold)
        for i in range(len(self.simData)):
            if i == num:
                redOrigin = tOrigin
            else:
                x = numpy.array(self.simData[i][0])
                y = numpy.array(self.simData[i][1])
                pylab.plot(x+tOrigin,y,colN)
                pylab.hold(True)
            # no need to increment tOrigin 
        x = numpy.array(self.simData[num][0])
        y = numpy.array(self.simData[num][1])
        pylab.plot(x+redOrigin,y,colY)

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
