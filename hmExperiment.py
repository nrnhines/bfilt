import hmm
import hmEnsemble
import numpy
import scipy.optimize

# A structure
def ch3bothdirs(tau01=2.,tau12=4.,nchannels=5):
    E = HME([])
    E.append(hmEnsemble.ch3Ensemble(V0=-65,V1=20,tau01=tau01,tau12=tau12,nchannels=nchannels))
    E.append(hmEnsemble.ch3Ensemble(V0=20,V1=-65,tau01=tau01,tau12=tau12,nchannels=nchannels))
    return E

# A structure
def ch3mix():
    E = HME([])
    E.append(hmEnsemble.ch3Ensemble(V0=-25.,V1=-20.,Vchar01=15.,Vchar12=15.))
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

def like4opt(param_values,param_names,known,structure,experiment):
    # like4opt(values=numpy.array([2,4]),
    # names=['tau01','tau12'],
    # known={'nchannels':5},
    # structure=ch3bothdirs)
    # experiment = HMEobject
    assert len(param_values)  == len(param_names)
    d = dict([(param_names[i],param_values[i]) for i in range(len(param_names))])
    d.update(known)
    S = structure(**d)
    L = S.likelihood(experiment)
    return -L

def like4eval(varied,fixed,structure,experiment):
    # for d in varied:
    #    fixed.update(d)
    S = structure(**fixed)
    L = S.likelihood(experiment)
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
            assert(isinstance(h,hmm.HMM))

    def sim(self,seeds):
        assert len(self) == len(seeds)
        for i in range(len(self)):
            self[i].sim(seeds[i])

    def simmed(self):
        s = True
        for M in self:
            s &= M.simmed
        return s

    def likelihood(self,simExperiment):
        assert simExperiment.simmed()
        total = 0
        for i in range(len(self)):
            total += self[i].likelihood(simExperiment[i].simData)
        return total

class fit(object):
    def __init__(self,guess,known,structure):
        self.guess = guess
        self.known = known
        self.structure = structure
        self.setted = False
        self.found = False

    def find(self,experiment):
        self.simExper = experiment
        vals = []
        names = []
        for key,v in guess.items():
            names.append(key)
            vals.append(v)
        values = numpy.array(vals)
        R = scipy.optimize.fmin_bfgs(like4opt,values,args=(names,known,structure,self.simExper),full_output=True,retall=True)
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

    def evalD(self,params,experiment=None,efun):
        if self.found:
            assert experiment == None  # self.simExper should already be defn
        else:
            self.simExper = experiment
        pnames = []
        ptuples = {EMPTY_SET}
        for k,v in guess.items():
            pnames.append(k)
            ptuples = ptuples {CARTESIAN PRODUCT} v
        L = []
        for p in ptuples:
            L = ptuples
            dictparam{FIGURE OUT dictionary(pnames,ptuples)}
            L = like4eval(**dictparam)
            zz.append(L)
        return (pnames, ptuples, zz)