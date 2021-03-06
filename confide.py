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
        self.funEval = self.profile2
        self.intFun = self.profile1
        self.convert = self.polar2x_y
        self.param0 = 0
        self.param1 = 1
        self.plotx = []
        self.ploty = []
        self.CI = []

    def setPlotParms(self,k0,k1):
        self.param0 = k0
        self.param1 = k1

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

    def meeker1(self,delta,k):
        CS = (delta**2)*self.H[k,k]
        return stats.chisqprob(CS,1) - self.alpha

    def profile1(self,delta,k):
        t = [k]
        n = []
        for i in range(int(self.MLE.size())):
            if not (i == k):
                n.append(i)
        A = self.H[numpy.ix_(t,t)]
        if len(n)>0:
            B = self.H[numpy.ix_(t,n)]
            D = self.H[numpy.ix_(n,n)]
            CS = (delta**2)*(A - B*D.I*B.T)[0,0]
        else:
            CS = (delta**2)*A[0,0]
        return stats.chisqprob(CS,2) - self.alpha

    def profile2(self,r,theta):
        k0 = self.param0
        k1 = self.param1
        x = r*math.cos(theta)
        y = r*math.sin(theta)
        delta = numpy.matrix([[x],[y]])
        t = [k0,k1]
        n = []
        for i in range(int(self.MLE.size())):
            if not (i ==k0 or i == k1):
                n.append(i)
        A = self.H[numpy.ix_(t,t)]
        if len(n)>0:
            B = self.H[numpy.ix_(t,n)]
            D = self.H[numpy.ix_(n,n)]
            CS = (delta.T*(A - B*D.I*B.T)*delta)[0,0]
        else:
            CS = (delta.T*A*delta)[0,0]
        return stats.chisqprob(CS,2) - self.alpha

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
        pe = self.intFun(d0,k)
        while pe > 0.0:
            d0 = d0*5
            pe = self.intFun(d0,k)
            assert(d0 < 5e7)
        delta = optimize.brentq(self.intFun,0,d0,(k,))
        return delta

    def getInterval(self,k):
        delta = self.find1(k)
        a = self.MLE[k] - delta
        b = self.MLE[k] + delta
        return (a,b)

    def CIs(self):
        self.CI = []
        for k in range(int(self.MLE.size())):
            (a,b) = self.getInterval(k)
            self.CI.append([a,b])
        return self.CI

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
            self.plotx.append(M[self.param0])
            self.ploty.append(M[self.param1])
        for rtheta in rt:
            self.plotx.append(M[self.param0] + rtheta[0]*math.cos(rtheta[1]))
            self.ploty.append(M[self.param1] + rtheta[0]*math.sin(rtheta[1]))
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

    def CIpoints(self, k0, k1):
        self.CIs()
        n = 50
        frac = (numpy.matrix(numpy.arange(0,n+1))/float(n))
        # print frac
        x0 = self.CI[k0][0] + frac*(self.CI[k0][1] - self.CI[k0][0])
        # x0 = numpy.arange(self.CI[k0][0],self.CI[k0][1],(self.CI[k0][1]-self.CI[k0][0])/50.0)
        y0 = []
        for x in x0.tolist():
            y0.append(self.MLE[k1])
        y1 = self.CI[k1][0] + frac*(self.CI[k1][1] - self.CI[k1][0])
        # y1 = numpy.arange(self.CI[k1][0],self.CI[k1][1],(self.CI[k1][1]-self.CI[k1][0])/50.0)
        x1 = []
        for y in y1.tolist():
            x1.append(self.MLE[k0])
        return (x0, y0, x1, y1)

    def plot(self,true=None):
        assert(len(self.plotx)>0)
        assert(len(self.plotx) == len(self.ploty))
        pyplot.ioff()
        pyplot.hold(False)
        pyplot.plot(self.plotx,self.ploty,'r*')
        pyplot.hold(True)
        MLE0 = self.MLE[self.param0]
        MLE1 = self.MLE[self.param1]
        # pyplot.plot([self.MLE[self.param0]],[self.MLE[self.param1]],'go')
        if not true == None:
            pyplot.plot([true[self.param0]],[true[self.param1]],'b+',markersize=20)
        # (x0, y0, x1, y1) = self.CIpoints(self.param0,self.param1)
        # print 'x0', x0, 'y0', y0, 'x1', x1, 'y1', y1
        # pyplot.plot(x0,y0,'go')
        # pyplot.plot(x1,y1,'go')
        L = self.CIs()
        delta0 = MLE0 - L[self.param0][0]
        delta1 = MLE1 - L[self.param1][0]
        pyplot.errorbar(MLE0,MLE1,delta1,delta0,'g')
        pyplot.title('Estimating Reaction Time Constants (3-State Channel)')
        pyplot.xlabel('Tau: Closed_1 To Open (msec)')
        pyplot.ylabel('Tau: Closed_2 To Closed_1 (msec)')
        pyplot.ion()
        pyplot.show()
