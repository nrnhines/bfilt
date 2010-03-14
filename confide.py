import numdifftools as nd
import math
from neuron import h
import scipy.stats as stats
import scipy.optimize as optimize
import matplotlib.pyplot as pyplot
import numpy

class Confide(object):
    def evalFun(self, p):
        saveParm = self.N.getParm()
        newParm = self.N.getParm()
        for i in range(len(p)):
            newParm.x[i] = p[i]
        self.N.setParm(newParm)
        L = self.N.likelihood()
        self.N.setParm(saveParm)
        return L

    def __init__(self,N,alpha=0.05):
        self.alpha = alpha
        self.N = N
        self.setMLE()
        self.funEval = self.hessEval
        self.convert = self.polar2x_y
        self.plotx = []
        self.ploty = []

    def Hessian(self):
        # Read current values of parameters (hopefully MLEs)
        parm = self.N.getParm()
        parmList = []
        for i in range(int(parm.size())):
            parmList.append(parm[i])
        # Create (likelihood) inline function
        LamFun = lambda p: self.evalFun(p)
        # Create Hessian (of likelihood) inline function
        HessFun = nd.Hessian(LamFun)
        # Evaluate Hessian and return
        return numpy.matrix(HessFun(parmList))

    def Gradient(self):
        # Read current values of parameters (hopefully MLEs)
        parm = self.N.getParm()
        parmList = []
        for i in range(int(parm.size())):
            parmList.append(parm[i])
        # Create (likelihood) inline function
        LamFun = lambda p: self.evalFun(p)
        # Create Hessian (of likelihood) inline function
        GradFun = nd.Gradient(LamFun)
        # Evaluate Hessian and return
        return numpy.matrix(GradFun(parmList))

    def setMLE(self):
        self.MLE = self.N.getParm()
        self.ML = self.N.likelihood()
        self.H = numpy.matrix(self.Hessian())
        self.InfoMatrix = self.H.I

    def return2MLE(self):
        self.N.setParm(self.MLE)

    def pValue(self,v):
        self.N.setParm(v)
        CS = 2.0*(self.N.likelihood() - self.ML)
        return stats.chisqprob(CS,v.size())

    def likeEval(self,r,theta):
        d = h.Vector(2)
        d.x[0] = self.MLE[0] + r*math.cos(theta)
        d.x[1] = self.MLE[1] + r*math.sin(theta)
        return self.pValue(d) - self.alpha

    def hessEval(self,r,theta):
        x = r*math.cos(theta)
        y = r*math.sin(theta)
        xy = numpy.matrix([[x],[y]])
        CS = xy.T*self.H*xy
        return stats.chisqprob(CS,2) - self.alpha

    def hess1(self,delta,k):
        CS = (delta**2)*self.H[k,k]
        return stats.chisqprob(CS,1) - self.alpha

    def hessTest(self,L):
        matrixList = L[:]
        for xy in matrixList:
            xy = numpy.matrix(xy).T
            xy[0] -= self.MLE[0]
            xy[1] -= self.MLE[1]
            CS = xy.T*self.H*xy
            print 'p-value:', stats.chisqprob(CS,2)

    def testxy(self, L):
        for xy in L:
            p = h.Vector(2)
            p.x[0] = xy[0]
            p.x[1] = xy[1]
            self.N.setParm(p)
            pv = pValue(nbf,p)
            print 'p-value', pv

    def polarFind(self,theta):
        self.return2MLE()
        r0 = 1
        pe = self.funEval(r0,theta)
        while pe > 0.0:
            r0 = r0*5
            pe = self.funEval(r0,theta)
            assert(r0 < 5e7)
        r = optimize.brentq(self.funEval,0,r0,(theta,))
        return r

    def find1(self,k):
        self.return2MLE()
        d0 = 1
        pe = self.hess1(d0,k)
        while pe > 0.0:
            d0 = d0*5
            pe = self.hess1(d0,k)
            assert(d0 < 5e7)
        delta = optimize.brentq(self.hess1,0,d0,(k,))
        return delta

    def getInterval(self,k):
        delta = self.find1(k)
        a = self.MLE[k] - delta
        b = self.MLE[k] + delta
        return (a,b)

    def polarInitial(self,r0,theta):
        self.return2MLE()
        z = self.funEval(r0,theta)
        if z > 0.0:
            a = r0
            b = r0*1.1
            while self.funEval(b,theta) > 0:
                b= b*1.1
                assert(b<1e9)
        elif z<0.0:
            b = r0
            a = r0*0.9
            while self.funEval(a,theta) < 0:
                a = a*0.9
                assert(a>1e-9)
        r = optimize.brentq(self.funEval,a,b,(theta,))
        return r

    def cont(self, theta0=0, points=125):
        r = self.polarFind(theta0)
        rt = [[r,theta0]]
        for i in range(points):
            theta = theta0 + 2*math.pi*(i+1)/points
            self.return2MLE()
            if self.funEval == self.likeEval:
                print 'point #', i
            r = self.polarInitial(r,theta)
            rt.append([r,theta])
        self.plotx = []
        self.ploty = []
        return self.convert(rt)

    def polar2x_y(self,rt,origin=None):
        self.plotx = []
        self.ploty = []
        if origin == None:
            M = self.MLE
        else:
            M = origin
            self.plotx.append(M[0])
            self.ploty.append(M[1])
        for rtheta in rt:
            self.plotx.append(M[0] + rtheta[0]*math.cos(rtheta[1]))
            self.ploty.append(M[1] + rtheta[0]*math.sin(rtheta[1]))
        return (self.plotx,self.ploty)

    def polar2polar(self,rt,origin=None):
        return rt

    def polar2xy(self,rt,origin=None):
        xy = []
        if origin == None:
            M = self.MLE
        else:
            M = origin
            xy.append([M[0],M[1]])
        for rtheta in rt:
            xy.append([M[0] + rtheta[0]*math.cos(rtheta[1]), M[1] + rtheta[0]*math.sin(rtheta[1])])
        return xy

    def plot(self,true=None):
        assert(len(self.plotx)>0)
        assert(len(self.plotx) == len(self.ploty))
        pyplot.ioff()
        pyplot.hold(False)
        pyplot.plot(self.plotx,self.ploty,'r*')
        pyplot.hold(True)
        pyplot.plot([self.MLE[0]],[self.MLE[1]],'go')
        if not true == None:
            pyplot.plot([true[0]],[true[1]],'b+')
        pyplot.ion()
        pyplot.show()