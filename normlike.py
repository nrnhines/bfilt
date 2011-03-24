import random
import math
import scipy.stats as stats

class normtcr(object):
    def __init__(self, seed=0, nsamp=100, mu=0, sig=1):
        self.seed = seed
        self.nsamp = nsamp
        self.mu = mu
        self.sig = sig
        self.R = random.Random()
        self.R.seed(self.seed)
        self.data = []
        for i in range(nsamp):
            self.data.append(self.R.normalvariate(self.mu,self.sig))
        self.mle = self.normmle()
        self.ml = self.like(self.mle)

    def normmle(self):
        muhat = 0
        for x in self.data:
            muhat += x
        muhat = muhat/len(self.data) #MLE = sample mean
        sigsquaredhat = 0
        for x in self.data:
            sigsquaredhat += (x-muhat)**2
        sigsquaredhat = sigsquaredhat/len(self.data)  #MLE not sample sig /n not /(n-1)
        sighat = math.sqrt(sigsquaredhat)
        return [muhat, sighat]

    def like(self,mutheta,sigma=None):
        if sigma==None:
            mu = mutheta[0]
            sig = mutheta[1]
        else:
            mu = mutheta
            sig = sigma
        mloglike = 0
        for x in self.data:
            mloglike += (1/2.)*math.log(2.*math.pi) + math.log(sig) + (x-mu)**2/(2.*sig**2)
        return mloglike

    def get_pValue(self,mutheta,sigma=None):
        minusPointLogLike = self.like(mutheta,sigma)
        minusMaxLogLike = self.ml
        size = 2
        CS = 2.0*(minusPointLogLike - minusMaxLogLike)
        pValue = stats.chisqprob(CS, size)
        return pValue

    def get_confidence(self,mutheta,sigma=None):
        return 100*(1-self.get_pValue(mutheta,sigma))

def evalpoint(x,y,efun):
    z = efun(x,y)
    return z

def eval2D(xx,yy,efun):
    zz = []
    for x in xx:
        zz.append([])
        for y in yy:
            zz[-1].append(efun(x,y))
            print 'x', x, 'y', y, 'z', zz[-1][-1]
    return zz

def save2D(fname,xx,yy,efun):
    zz = eval2D(xx,yy,efun)
    f = open(fname+"_x.txt","w")
    for x in xx:
        f.write(str(x)+' ')
    f.close()
    f = open(fname+"_y.txt","w")
    for y in yy:
        f.write(str(y)+' ')
    f.close()
    f = open(fname+"_z.txt","w")
    for zr in zz:
        for z in zr:
            f.write(str(z)+' ')
    f.close()

