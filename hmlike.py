import hmm
import hmEnsemble

class HML(object):
    def __init__(self,tau01=2.,tau12=4.,N=5):
        self.tau01_true = tau01
        self.tau12_true = tau12
        self.N_true = N
        self.M_true = hmEnsemble.ch3Ensemble(tau01=tau01,tau12=tau12,nchannels=N)

    def sim(self,seed=0):
        self.seed = seed
        self.M_true.sim(seed=self.seed)

    def ch3like(self,tau01,tau12,N):
        M_assumed = hmEnsemble.ch3Ensemble(tau01=tau01,tau12=tau12,nchannels=N)
        L = M_assumed.likelihood(self.M_true.simData)
        return L

    def evalsave2D(fname,xx,yy,N):
        self.eval2D(xx,yy,N)
        self.save2D(fname)

    def eval2D(self,xx,yy,N_assumed):
        self.xx = xx
        self.yy = yy
        self.N_assumed = N_assumed
        self.zz = []
        for x in self.xx:
            self.zz.append([])
            for y in self.yy:
                self.zz[-1].append(self.ch3like(x,y,self.N_assumed))
                print 'x', x, 'y', y, 'z', self.zz[-1][-1]

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
