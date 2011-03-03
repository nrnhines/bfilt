import random
import math
import pylab
import numpy
import scipy.stats as stats

class normtcr(object):
    def __init__(self, seed=0, nsamp=100, mu=0, sig=1):
        save = True # True until all recursive computations defined
        self.seed = seed
        self.nsamp = nsamp
        self.mu = mu
        self.sig = sig
        self.R = random.Random()
        self.R.seed(self.seed)
        self.data = self.gen(save)
        self.calc()

    def calc(self):
        self.mle = self.normmle()
        self.ml = self.like(self.mle)
        self.pValue = self.get_pValue(self.mu,self.sig)
        self.mubias = self.mle[0] - self.mu
        self.sigbias = self.mle[1] - self.sig
        self.normbias = math.sqrt(self.mubias**2 + self.sigbias**2)

    def reseed(self,seed):
        save = True
        self.seed = seed
        self.R.seed(self.seed)
        self.data = self.gen(save)
        self.calc()

    def gen(self,save):
        if save:
            data = []
        else:
            data = None
        for i in range(self.nsamp):
            d = self.R.normalvariate(self.mu,self.sig)
            # Commented out: recursively compute muhat
            #~ if i == 0:
                #~ self.muhat = d
            #~ else:
                #~ self.muhat = (i/(i+1.))*self.muhat + (1./(i+1.))*d
            #~ if save:
                #~ data.append(d)
            data.append(d)
        return data

    def normmle(self):
        muhat = 0
        for x in self.data:
            muhat += x
        muhat = muhat/len(self.data) #MLE = sample mean
        sigsquaredhat = 0
        for x in self.data:
            sigsquaredhat += (x-muhat)**2
        sigsquaredhat = sigsquaredhat/len(self.data)  #MLE is not sample sig: /n not /(n-1)
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

class cexperiment(object):
    def __init__(self, seed0=0, nruns=10, reps=1000, nsamp=4, mu=0, sig=1):
        self.N = normtcr(seed=seed0,nsamp=nsamp,mu=mu,sig=sig)
        self.seed0=seed0
        self.nruns=nruns
        self.reps = reps
        self.nsamp=nsamp
        self.mu=mu
        self.sig=sig
        self.covers = []
        self.cutoff = []
        for i in range(nruns):
            s = seed0+i
            self.N.reseed(s)
            self.NE = nexperiment(seed0=s*10000,nruns=self.reps,nsamp=[self.nsamp],mu=self.N.mle[0],sig=self.N.mle[1])
            self.covers.append(self.NE.covers[0])
            self.NE.pvs[0].sort()
            self.cutoff.append(self.NE.pvs[0][reps/20])

class nexperiment(object):
    def __init__(self, seed0=0, nruns=100, nsamp=[4,8,16,32,64,128,256,512], mu=0, sig=1):
        N = normtcr(seed=seed0,nsamp=nsamp[0],mu=mu,sig=sig)
        self.nsamp = nsamp
        self.nruns = nruns
        self.mu = mu
        self.sig = sig
        self.seed0 = seed0
        self.mbias = []
        self.sembias = []
        self.pvs = []
        self.alpha = 0.05
        for s in nsamp:
            N.nsamp=s
            normbias = []
            self.pvs.append([])
            for r in range(nruns):
                N.reseed(seed0+r)
                normbias.append(N.normbias)
                self.pvs[-1].append(N.pValue)
            seed0 += nruns
            (mb,semb) = self.meansem(normbias)
            self.mbias.append(mb)
            self.sembias.append(semb)
        self.covers = self.newcovers(self.alpha)

    def newcovers(self,alpha):
        covers = []
        for pv in self.pvs:
            c = 0
            for p in pv:
                if p>self.alpha:
                    c += 1
            covers.append(c)
        return covers

    def meansem(self,L):
        sumL = 0.0
        for l in L:
            sumL += l
        m = sumL/len(L)
        sumS = 0.0
        for l in L:
            sumS += (l - m)**2
        s = math.sqrt(sumS/(len(L)-1))
        sem = s/math.sqrt(len(L))
        return (m,sem)

    def pvhist(self,ns=0,nbins=20):
        print "nsamp", self.nsamp[ns]
        pylab.ion()
        pylab.hist(numpy.array(self.pvs[ns]),nbins)

    def pvdump(self,fname):
        f = open(fname,'w')
        for j in range(len(self.nsamp)):
            f.write("samp"+str(self.nsamp[j])+' ')
        f.write('\n')
        for j in range(len(self.pvs[0])):
            for k in range(len(self.pvs)):
                f.write(str(self.pvs[k][j])+' ')
            f.write("\n")
        f.close()

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

