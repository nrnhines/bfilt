import numpy
import random

class Regress(object):
    def __init__(self,x,y):
        self.Rand = random.Random()
        self.Rand.seed(0)
        self.px = x
        self.py = y
        self.dx = []
        self.dy = []
        self.sigma = None
        self.seed = 0
        self.n = 0
        assert(len(x) == len(y))
        self.trueOrder = len(x)-1
        self.order = self.trueOrder
        a = numpy.polyfit(x,y,self.trueOrder)
        self.poly = []
        for i in range(len(a)):
            self.poly.append(a[-i-1])

    def badpolyeval(self,x):
        y = 0
        for i in range(len(self.poly)):
            y += self.poly[i]*(x**i)
        return y

    def polyeval(self, x):
        y = 0
        rmod = self.poly[:]
        rmod.reverse()
        for coef in rmod:
            y = y * x + coef
        return y

    def genData(self,n=50,sigma=0.01,seed=0):
        self.dx = []
        self.dy = []
        self.n = n
        self.sigma = sigma
        self.seed = seed
        self.Rand.seed(seed)
        for k in range(n):
            self.dx.append(self.Rand.random())
            self.dy.append(self.Rand.normalvariate(0,sigma))
            self.dy[-1] += self.polyeval(self.dx[-1])

    def dataPairs(self):
        assert(len(self.dx) == len(self.dy))
        xy = []
        for i in range(len(self.dx)):
            xy.append([self.dx[i], self.dy[i]])
        return xy

    def polyPairs(self,t0=0.0,t1=1.0,dt=0.01):
        mxy = []
        for t in numpy.arange(t0,t1,dt):
            mxy.append([t, self.polyeval(t)])
        return mxy

    def likelihood():
        loglike = NumPoints*log(1/(pdata.NoiseLevel*sqrt(2*pi))) - 1/(2*pdata.NoiseLevel^2)*S.normr^2;

    def getParm():

    def setParm():

    def findMLE():
        a = numpy.polyfit(self.dx,self.dy,self.order)
        self.mod = []
        for i in range(len(a)):
            self.mod.append(a[-i-1])

    def setOrder(order=None):
        if order == None:
            self.order = self.trueOrder
        else
            self.order = order
