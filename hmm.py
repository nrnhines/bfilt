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

def ch3Qv(V, tau01, tau12):
    inf01 = 1./(1. + math.exp(1.*(-20. - V)))
    inf12 = 1./(1. + math.exp(1.*(-25. - V)))
    alpha01 = tau01*inf01
    beta01 = tau01-alpha01
    alpha12 = tau12*inf12
    beta12 = tau12-alpha12
    Q = ch3Q(alpha01, beta01, alpha12, beta12)
    return Q

def ch3hmm(V0=-65, V1=20, tau01=2, tau12=4,sigma=0.001):
    Q0 = ch3Qv(V0, tau01, tau12)
    pstates = equilibrium(Q0)
    output = [0.0, 0.0, 1.0]
    Q = ch3Qv(V1, tau01, tau12)
    H = HMM(pstates,output,Q,sigma)
    return H

def equilibrium(Q):
    (V,D) = numpy.linalg.eig(Q.T)
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
    tol = 1e-6
    print 'pstates', pstates
    for i in range(len(pstates)):
        assert(pstates[i] >= 0.0)
        assert(pstates[i] <= 1.0)
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

    def sim(self, seed=0, dt=0.1, tstop=20):
        self.seed = seed
        self.dt = dt
        self.tstop = tstop
        self.trans = scipy.linalg.expm(dt*self.Q)
        tol = 1e-7
        for i in range(self.nstates):
            rowsum = 0
            for j in range(self.nstates):
                rowsum += self.trans[i,j]
            assert math.fabs(rowsum - 1.0) < tol
        assert(tstop>0.0)
        print "Checks out..."
        self.R.seed(seed)
        self.nsamples = int(math.ceil(tstop/dt))
        self.simStates = []
        self.simStates.append(self.select(self.init,0))
        for i in range(self.nsamples-1):
            self.simStates.append(self.select(self.trans,self.simStates[-1]))
        # print 'States as a function of discrete time:', self.simStates
        self.simOut = []
        for s in self.simStates:
            self.simOut.append(self.output[s])
        # print 'Output with out noise:', self.simOut
        self.simData = []
        for o in self.simOut:
            self.simData.append(o + self.R.normalvariate(0,self.sigma))
        # print 'simData:', self.simData
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

    def predict(self,inter):
        # print 'time', self.time[-1], 'pmf', inter
        self.saveErrorBars(inter,self.time[-1])
        for i in range(self.skip):
            inter = inter*self.trans
            self.saveErrorBars(inter,self.time[-1]+self.dt)
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

    def likelihood(self, fitData, skip):
        self.fitData = fitData
        self.skip = skip
        self.transfit = scipy.linalg.expm(self.dt*self.skip*self.Q)
        self.initializeErrorBars()
        pmf = self.init
        sll = 0
        for i in range(self.skip,len(self.fitData),self.skip):
            # print i, i*self.dt #, self.fitData[i]
            pre = self.predict(pmf)
            (pmf, lk) = self.update(self.fitData[i], pre)
            sll += math.log(lk)
        self.liked = True
        return sll

    def plot(self):
        assert(self.liked)
        x = numpy.arange(0,self.dt*self.nsamples,self.dt)
        y = numpy.array(self.fitData)
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
