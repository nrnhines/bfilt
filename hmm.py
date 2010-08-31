import pylab
import numpy
import random
import math
import scipy
import scipy.linalg

class HMM(object):
    def __init__(self, param):
        dt = 1.0
        self.dt = dt
        skip = 20
        self.skip = skip
        init = [1.0,0.0,0.0]
        output = [1.0,0.0,0.0]
        sigma = 0.01
        self.nstates = len(init)
        self.Q = self.makeQ(param[0], param[1], param[2], param[3])
        self.trans = scipy.linalg.expm(dt*self.Q)
        self.transfit = scipy.linalg.expm(dt*skip*self.Q)
        self.init = numpy.matrix(init)
        self.output = output
        self.sigma = sigma
        self.simStates = None
        self.simOut = None
        self.nsamples = None
        self.Data = None
        self.center = None
        self.width = None
        self.R = random.Random()
        tol = 1e-7
        assert self.init.shape[1] == self.trans.shape[0]
        assert self.init.shape[1] == self.trans.shape[1]
        assert self.init.shape[0] == 1
        assert len(output) == self.init.shape[1]
        assert math.fabs(sum(init) - 1.0) < tol
        for i in range(self.nstates):
            rowsum = 0
            for j in range(self.nstates):
                rowsum += self.trans[i,j]
            assert math.fabs(rowsum - 1.0) < tol
        print "Checks out..."

    def makeQ(self, alpha, beta, gamma, delta):
        Q = numpy.matrix([[-alpha, alpha, 0], [beta, -(beta + gamma), gamma], [0, delta, -delta]])
        return Q

    def select(self,mat,row):
        p = self.R.random()
        rowsum = 0
        for i in range(mat.shape[1]):
            rowsum += mat[row,i]
            if p < rowsum:
                return i
        assert False

    def sim(self, seed, nsamples):
        self.R.seed(seed)
        self.nsamples = nsamples
        self.simStates = []
        self.simStates.append(self.select(self.init,0))
        for i in range(nsamples-1):
            self.simStates.append(self.select(self.trans,self.simStates[-1]))
        # print 'States as a function of discrete time:', self.simStates
        self.simOut = []
        for s in self.simStates:
            self.simOut.append(self.output[s])
        # print 'Output with out noise:', self.simOut
        self.Data = []
        for o in self.simOut:
            self.Data.append(o + self.R.normalvariate(0,self.sigma))
        # print 'Data:', self.Data

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

    def predict(self,start):
        inter = start
        self.saveErrorBars(inter,self.time[-1])
        for i in range(self.skip):
            inter = inter*self.trans
            self.saveErrorBars(inter,self.time[-1]+self.dt)
        return inter

    def collect(self,data):
        col = []
        for o in self.output:
            col.append(self.normpdf(data - o,0,self.sigma))
        return col

    def weigh(self,weight,pmf):
        new = []
        for i in range(len(weight)):
            new.append(weight[i]*pmf[0,i])
        return new

    def marginal(self,new):
        return(sum(new))

    def update(self,new,marg):
        return numpy.matrix(new)/marg

    def likelihood(self):
        self.initializeErrorBars()
        pmf = self.init
        sll = 0
        for i in range(0,len(self.Data),self.skip):
            d = self.Data[i]
            pre = self.predict(pmf)
            col = self.collect(d)
            new = self.weigh(col, pre)
            lk = self.marginal(new)
            pmf = self.update(new,lk)
            sll += math.log(lk)
        return sll

    def plot(self):
        x = numpy.arange(0,self.dt*self.nsamples,self.dt)
        y = numpy.array(self.Data)
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
