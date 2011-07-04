import hmm
import hmEnsemble
import numpy
import scipy.optimize
import scipy.stats as stats
import numdifftools as nd
import copy

def ch3like4opt(param,V0,V1,N,Data):
    tau01 = abs(param[0])
    tau12 = abs(param[1])
    M_assumed = hmEnsemble.ch3Ensemble(V0=V0,V1=V1,tau01=tau01,tau12=tau12,nchannels=N)
    try:
        L = M_assumed.likelihood(Data)
    except:
        print "Out of range: tau01",p[0], "tau12",p[1]
        L = numpy.nan
    return -L

class HMMLike(object):
    def __init__(self,System,Model):
        self.System = System
        self.Model = Model
    
    def sim(self, seeds=[0], dt=0.1, tstops=[20]):
        System.sim(seeds,dt,tstops)
    
    
            
class HML(object):
    def __init__(self,V0=-65,V1=20,tau01=2.,tau12=4.,N=5):
        self.V0 = V0
        self.V1 = V1
        self.tau01_true = tau01
        self.tau12_true = tau12
        self.N_true = N
        self.M_true = hmEnsemble.ch3Ensemble(V0=V0,V1=V1,tau01=tau01,tau12=tau12,nchannels=N)
        self.ml = None

    def sim(self,seeds=[0]):
        self.seeds = seeds
        self.M_true.sim(seeds=self.seeds)

    def ch3like(self,tau01,tau12,N):
        M_assumed = hmEnsemble.ch3Ensemble(V0=self.V0,V1=self.V1,tau01=tau01,tau12=tau12,nchannels=N)
        #try:
        L = M_assumed.likelihood(self.M_true.simData)
        #except:
            #print "Out of range: tau01", tau01, "tau12", tau12
            #L = numpy.nan
        return L

    # def setN(self,N):
    #      self.N_assumed = N

    def find(self,tau01,tau12,N):
        p0 = numpy.array([tau01,tau12])
        # xopt,fopt,gopt,Bopt,funcalls,gradcalls,warnflag,allvecs =
        return scipy.optimize.fmin_bfgs(ch3like4opt,p0,args=(self.V0,self.V1,N,self.M_true.simData),full_output=True,retall=True)

    def findmle(self,tau01,tau12,N):
        R = self.find(tau01,tau12,N)
        self.mle = R[0]
        self.ml = -R[1]
        self.fopt = R[1]
        self.gopt = R[2]
        self.Bopt = R[3]
        self.funcalls = R[4]
        self.gradcalls = R[5]
        self.warnflag = R[6]
        self.allvecs = R[7]
        return self.mle

    def p_true(self):
        self.truelike = self.M_true.likelihood(self.M_true.simData)
        self.findmle(self.tau01_true,self.tau12_true,self.N_true)
        self.CS = 2.0*(self.ml - self.truelike)  #plus log-likelihood
        numparam = 2
        self.pValue = stats.chisqprob(self.CS,numparam)
        return self.pValue

    def ch3p(self,tau01,tau12,N_assumed):
        # Must compute ml first
        like = self.ch3like(tau01,tau12,N_assumed)
        numparam = 2
        return stats.chisqprob(2.*(self.ml - like),numparam) # plus log-likelihood

    def ch3c(self,tau01,tau12,N_assumed):
        p = self.ch3p(tau01,tau12,N_assumed)
        return 100.*(1.-p)

    def evallike(self,param,N_assumed):
        return self.ch3like(param[0],param[1],N_assumed)

    def Hessian(self,tau01,tau12,N):
        # Create (likelihood) inline function
        LamFun = lambda p: self.evallike(p,N)
        # Create Hessian (of likelihood) inline function
        HessFun = nd.Hessian(LamFun)
        # Evaluate Hessian and return
        return numpy.matrix(HessFun([tau01,tau12]))

    def evalsave2D(fname,xx,yy,N):
        self.eval2D(xx,yy,N)
        self.save2D(fname)

    def eval2D(self,xx,yy,N_assumed,efun=None):
        if efun == None:
            efun=self.ch3like
        self.xx = xx
        self.yy = yy
        self.N_assumed = N_assumed
        self.zz = []
        for x in self.xx:
            self.zz.append([])
            for y in self.yy:
                self.zz[-1].append(efun(x,y,self.N_assumed))
                print 'x', x, 'y', y, 'z', self.zz[-1][-1]

    def eval2Dp(self,xx,yy,N_assumed):
        self.findmle(self.tau01_true,self.tau12_true,self.N_true)
        self.eval2D(xx,yy,N_assumed,self.ch3p)

    def eval2Dc(self,xx,yy,N_assumed):
        self.findmle(self.tau01_true,self.tau12_true,self.N_true)
        self.eval2D(xx,yy,N_assumed,self.ch3c)

    def save2D(self,fname):
        self.fname = fname
        f = open(self.fname+"_x.txt","w")
        for x in self.xx:
            f.write(str(x)+' ')
        f.close()
        f = open(self.fname+"_y.txt","w")
        for y in self.yy:
            f.write(str(y)+' ')
        f.close()
        f = open(self.fname+"_z.txt","w")
        for zr in self.zz:
            for z in zr:
                f.write(str(z)+' ')
        f.close()
