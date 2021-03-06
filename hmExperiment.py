import hmm
import hmEnsemble
import numpy
import pylab
import scipy.optimize
import scipy.stats
import copy

# A structure
def ch3up(tau01=2.,tau12=4.,nchannels=5,Vchar01=1.,Vchar12=1.):
    E = HME([])
    E.append(hmEnsemble.ch3Ensemble(V0=-65,V1=20,tau01=tau01,tau12=tau12,Vchar01=Vchar01,Vchar12=Vchar12,nchannels=nchannels))
    return E

# A structure
def ch3bothdirs(tau01=2.,tau12=4.,nchannels=5):
    E = HME([])
    E.append(hmEnsemble.ch3Ensemble(V0=-65,V1=20,tau01=tau01,tau12=tau12,nchannels=nchannels))
    E.append(hmEnsemble.ch3Ensemble(V0=20,V1=-65,tau01=tau01,tau12=tau12,nchannels=nchannels))
    return E

# A structure
def ch3mix():
    E = HME([])
    E.append(hmEnsemble.ch3Ensemble(tau01=tau01,tau12=tau12,V0=-25.,V1=-20.,Vchar01=15.,Vchar12=15.,nchannels=nchannels))
    return E

# A structure
def ch3chn(tau01=2.,tau12=4.,nchannels=5,Vchar01=1.,Vchar12=1.):
    E = HME([])
    E.append(hmEnsemble.ch3EnsemChain(tau01=tau01,tau12=tau12,nchannels=nchannels,Vchar01=Vchar01,Vchar12=Vchar12))
    return E

# A structure
def hybridprotocol(tau01=2.,tau12=4.,nchannels=5):
    E = HME([])
    E.append(hmEnsemble.ch3EnsemChain(tau01=tau01,tau12=tau12,nchannels=nchannels))
    E.append(hmEnsemble.ch3Ensemble(tau01=tau01,tau12=tau12,V0=-25.,V1=-20.,Vchar01=15.,Vchar12=15.,nchannels=nchannels))
    return E

def testfind():
    E_true = ch3bothdirs()
    E_true.sim([[1,2],[4,5]])
    guess = dict(tau01=2.,tau12=4.)
    known = dict(nchannels=5)
    structure = ch3bothdirs
    data = E_true
    MLE = E_true.find(guess, known, structure)
    print "MLE", MLE

def like4opt(param_values,param_names,known,structure,SysWData):
    # like4opt(
    #   values=numpy.array([2,4]),
    #   names=['tau01','tau12'],
    #   known={'nchannels':5},
    #   structure=ch3bothdirs)
    #   sysWData=HMEobject
    # )
    assert len(param_values) == len(param_names)
    d = dict([(param_names[i],param_values[i]) for i in range(len(param_names))])
    d.update(known)
    S = structure(**d)
    L = S.likelihood(SysWData)
    return -L

def like4eval(params,structure,SysWData):
    S = structure(**params)
    L = S.likelihood(SysWData)
    return L

class HME(object):
    def __init__(self,protocol):
        self.integrity(protocol)
        self.protocol = protocol

    def __len__(self):
        return len(self.protocol)

    def __getitem__(self,key):
        return self.protocol[key]

    def __setitem__(self,key,value):
        self.integrity([value])
        self.protocol[key] = value

    def __delitem__(self,key):
        del self.protocol[key]

    def __iter__(self):
        return self.protocol.__iter__()

    def __add__(self,other):
        return HME(self.protocol + other.protocol)

    def append(self,other):
        self.integrity([other])
        self.protocol.append(other)

    def integrity(self,list=None):
        if list == None:
            list = self.protocol
        for h in list:
            assert(isinstance(h,hmm.HMM) or isinstance(h,hmm.HMMChain))

    def sim(self,seeds,dt,tstops):
        if tstops == None:
            tstops = [None]*len(self)
        if not type(tstops) is list:
            tstops = [tstops]*len(self)
        assert len(self) == len(seeds)
        for i in range(len(self)):  # i ranges over the experiments
            self[i].sim(seeds[i],dt,tstops[i])

    def simmed(self):
        s = True
        for M in self:
            s &= M.simmed
        return s
    
    def simplot(self,num):
        assert self.simmed()
        fig = 0
        for E in self.protocol:
            pylab.figure(fig)
            E.simplot(num)
            fig += 1

    def likelihood(self,simExperiment):
        assert simExperiment.simmed
        total = 0
        for i in range(len(self)):
            total += self[i].likelihood(simExperiment[i])
        return total

class fit(object):
    def __init__(self,structure,guess,known):
        self.guess = guess
        self.known = known
        self.structure = structure
        self.found = False
        self.simmed = False

    def like(self,point,SysWData):
        p = dict()
        p.update(self.known)
        p.update(self.guess)
        p.update(point)
        S = self.structure(**p)
        L = S.likelihood(SysWData)
        return L
    
    def sim(self,system,true,seeds=[0],dt=.1,tstops=None):
        self.system = system
        self.true = true
        self.SysWData = system(**true)
        self.SysWData.sim(seeds,dt,tstops)
        self.simmed = True

    def find(self):
        assert(self.simmed)
        vals = []
        names = []
        for key,v in self.guess.items():
            names.append(key)
            vals.append(v)
        values = numpy.array(vals)
        R = scipy.optimize.fmin_bfgs(like4opt,values,args=(names,self.known,self.structure,self.SysWData),full_output=True,retall=True)
        self.MLE = R[0]
        self.ML = -R[1]
        self.fopt = R[1]
        self.gopt = R[2]
        self.Bopt = R[3]
        self.funcalls = R[4]
        self.gradcalls = R[5]
        self.warnflag = R[6]
        self.allvecs = R[7]
        self.found = True
        return self.MLE

    def simplot(self,num=0):
        assert(self.simmed)
        self.SysWData.simplot(num)
    
    # The following routine is not functional
    def evalD(self,params,efun=like4eval):
        pnames = params.keys()
        pfactorlist = params.values()
        zz = {}
        for ptuples in myutil.cartesian_product(*pfactorlist):  # myutil not written
            paramset = dict([(pnames[i],ptuples[i]) for i in range(len(ptuples))])
            L = efun(paramset,self.structure,self.SysWData)
            zz.update({ptuple:L})
        return (pnames, zz)

    def pValue(self,z,ml=None,nparam=None):
        if ml == None:
            ml = self.ML
        if nparam == None:
            nparam = len(self.guess)
        CS = 2.0*(ml-z)  #plus log-likelihood
        p = scipy.stats.chisqprob(CS,nparam)
        return p
    
    def cValue(self,z,ml=None,nparam=None):
        p = self.pValue(z,ml,nparam)
        return 100*(1-p)
    
    def save4plot(self,fname,params,xname="tau01",yname="tau12",base={},efun=like4eval):
        mle = self.find()
        self.fname = fname
        f = open(self.fname+"_x.txt","w")
        for x in params[xname]:
            f.write(str(x)+' ')
        f.close()
        f = open(self.fname+"_y.txt","w")
        for y in params[yname]:
            f.write(str(y)+' ')
        f.close()
        f = open(self.fname+"_z.txt","w")
        g = open(self.fname+"_c.txt","w")
        p = copy.deepcopy(params)
        p.update(base)
        for x in params[xname]:
            p.update({xname:x})
            for y in params[yname]:
                print xname, x, yname, y
                p.update({yname:y})
                # Changed plans from: for pn in pnzz[0], pnzz = (pnames,zz)
                z = efun(p,self.structure,self.SysWData)
                f.write(str(z)+' ')
                c = self.cValue(z)
                g.write(str(c)+' ')
        f.close()
        g.close()
