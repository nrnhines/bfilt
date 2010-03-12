import numpy
import random
import math
from neuron import h

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

    def eval(self, p, x):
        y = 0
        rp = p[:]
        rp.reverse()
        for coef in rp:
            y = y * x + coef
        return y

    def polyeval(self, x):
        return self.eval(self.poly,x)

    def modeval(self, x):
        return self.eval(self.mod,x)

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

    def modPairs(self,t0=0.0,t1=1.0,dt=0.01):
        mxy = []
        for t in numpy.arange(t0,t1,dt):
            mxy.append([t, self.modeval(t)])
        return mxy

    def likelihood(self):
        RSS = 0
        assert(len(self.dx) == len(self.dy))
        for i in range(len(self.dx)):
            RSS += (self.dy[i] - self.modeval(self.dx[i]))**2
        loglike = self.n*math.log(1.0/(self.sigma*math.sqrt(2.0*math.pi))) - 1.0/(2.0*self.sigma**2)*RSS;
        return -loglike

    def getParm(self):
        Parm = h.Vector(len(self.mod))
        for i in range(len(self.mod)):
            Parm.x[i] = self.mod[i]
        return Parm

    def setParm(self,v):
        self.order = int(v.size()) - 1
        for i in range(len(self.mod)):
            self.mod[i] = v[i]

    def getMod(self):
        return self.mod

    def setMod(self,m):
        self.order = len(m) - 1
        self.mod = m

    def findMLE(self,order=None):
        if order == None:
            order = self.order
        a = numpy.polyfit(self.dx,self.dy,order)
        self.mod = []
        for i in range(len(a)):
            self.mod.append(a[-i-1])

    def findTrue(self):
        self.mod = self.poly

    def setOrder(order=None):
        if order == None:
            self.order = self.trueOrder
        else:
            self.order = order
