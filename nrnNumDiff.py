import numdifftools as nd
import math
from neuron import h
import scipy.stats as stats
import scipy.optimize as optimize

def evalFun(nbf, p):
    saveParm = nbf.getParm()
    newParm = nbf.getParm()
    for i in range(len(p)):
        newParm.x[i] = p[i]
    nbf.setParm(newParm)
    L = nbf.likelihood()
    nbf.setParm(saveParm)
    return L

def Hessian(nbf):
    # Read current values of parameters (hopefully MLEs)
    parm = nbf.getParm()
    parmList = []
    for i in range(int(parm.size())):
        parmList.append(parm[i])
    # Create (likelihood) inline function
    LamFun = lambda p: evalFun(nbf,p)
    # Create Hessian (of likelihood) inline function
    HessFun = nd.Hessian(LamFun)
    # Evaluate Hessian and return
    return HessFun(parmList)

def Gradient(nbf):
    # Read current values of parameters (hopefully MLEs)
    parm = nbf.getParm()
    parmList = []
    for i in range(int(parm.size())):
        parmList.append(parm[i])
    # Create (likelihood) inline function
    LamFun = lambda p: evalFun(nbf,p)
    # Create Hessian (of likelihood) inline function
    GradFun = nd.Gradient(LamFun)
    # Evaluate Hessian and return
    return GradFun(parmList)

def LikePerturbed(nbf,perturbList,delta):
    saveParm = nbf.getParm()
    newParm = nbf.getParm()
    for i in range(len(perturbList)):
        newParm.x[i] = newParm[i] + delta*perturbList[i]
    nbf.setParm(newParm)
    L = nbf.likelihood()
    nbf.setParm(saveParm)
    return L

def test4(Parm):
    LamFun = lambda p: p**4
    HessFun = nd.Hessian(LamFun)
    return HessFun(Parm)

def rememberMLE(nbf):
    global ML, MLE
    MLE = nbf.getParm()
    ML = nbf.likelihood()

def pMinusAlpha(nbf,v,alpha):
    global ML
    # print 'before',v
    # print 'here', v[0]
    # print 'there', v[1]
    nbf.setParm(v)
    CS = 2.0*(nbf.likelihood() - ML)
    return stats.chisqprob(CS,v.size()) - alpha

def perturbedEval(pert,nbf,dir,alpha):
    global MLE
    # print 'One', MLE
    # print 'MLE[0]', MLE[0]
    # print 'MLE[1]', [1]
    v = MLE.c()
    d = dir.c()
    # print 'Two', v
    # print '2:v[0]', v[0]
    # print '2:v[1]', v[1]
    d.mul(pert)
    # print 'dir.mul', dir
    v.add(d)
    # print 'Four', v
    # print '4:here', v[0]
    # print '4:there', v[1]
    return pMinusAlpha(nbf,v,alpha)

def findSignificant(nbf,alpha):
    global ML, MLE
    nbf.setParm(MLE)
    dir = h.Vector(2)
    dir.x[0] = 1
    dir.x[1] = 0
    pert = 1
    # print 'val', perturbedEval(pert,nbf,dir,alpha)
    # print 'truth', perturbedEval(pert,nbf,dir,alpha) > 0.0
    while perturbedEval(pert,nbf,dir,alpha) > 0.0:
        pert = pert*5
        print 'pert* = ', pert
        assert(pert < 500000)
    # print 'f(b)', perturbedEval(pert,nbf,dir,alpha)
    # print 'f(a)', perturbedEval(0,nbf,dir,alpha)
    # print 'pert', pert
    p = optimize.brentq(perturbedEval,0,pert,(nbf,dir,alpha))
    dir.mul(p)
    return dir


