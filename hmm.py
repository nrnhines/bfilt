import numpy
import random
import math

class HMM(object):
    def __init__(self, trans, init, output, sigma):
        self.nstates = len(init)
        self.trans = trans
        self.init = numpy.matrix(init)
        self.output = output
        self.sigma = sigma
        self.simStates = None
        self.simOut = None
        self.Data = None
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
        self.simStates = []
        self.simStates.append(self.select(self.init,0))
        for i in range(nsamples-1):
            self.simStates.append(self.select(self.trans,self.simStates[-1]))
        print 'States as a function of discrete time:', self.simStates
        self.simOut = []
        for s in self.simStates:
            self.simOut.append(self.output[s])
        print 'Output with out noise:', self.simOut
        self.Data = []
        for o in self.simOut:
            self.Data.append(o + self.R.normalvariate(0,self.sigma))
        print 'Data:', self.Data

    def normpdf(self,x,m,sigma):
        return (1/(math.sqrt(2*math.pi)*sigma))*math.exp(-(x-m)*(x-m)/(2*sigma*sigma))

    def predict(self,start):
        return start*self.trans

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
        pmf = self.init
        sll = 0
        for d in self.Data:
            pre = self.predict(pmf)
            col = self.collect(d)
            new = self.weigh(col, pre)
            lk = self.marginal(new)
            pmf = self.update(new,lk)
            sll += math.log(lk)
        return sll