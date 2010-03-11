import numdifftools as nd
import math
from neuron import h
import scipy.stats as stats
import scipy.optimize as optimize
import numpy

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
    global H, ML, MLE
    MLE = nbf.getParm()
    ML = nbf.likelihood()
    H = Hessian(nbf)

def return2MLE(nbf):
    global MLE
    nbf.setParm(MLE)

def pValue(nbf,v):
    global ML
    nbf.setParm(v)
    CS = 2.0*(nbf.likelihood() - ML)
    return stats.chisqprob(CS,v.size())

def pMinusAlpha(nbf,v,alpha):
    global ML
    # print 'before',v
    # print 'here', v[0]
    # print 'there', v[1]
    # print 'likelihood', nbf.likelihood()
    # print 'ML', ML
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

def likeEval(r,theta,nbf,alpha):
    global MLE
    d = h.Vector(2)
    d.x[0] = MLE[0] + r*math.cos(theta)
    d.x[1] = MLE[1] + r*math.sin(theta)
    return pMinusAlpha(nbf,d,alpha)

def hessEval(r,theta, nbf, alpha=0.0):
    global H, MLE
    x = r*math.cos(theta)
    y = r*math.sin(theta)
    HI = numpy.matrix(H).I
    xy = numpy.matrix([[x],[y]])
    CS = xy.T*HI*xy
    return stats.chisqprob(CS,2)-alpha

def hessTest(List):
    global H,ML
    HI = numpy.matrix(H).I
    copyList = List[:]
    for xy in copyList:
        xy = numpy.matrix(xy).T
        xy[0] -= MLE[0]
        xy[1] -= MLE[1]
        CS = xy.T*HI*xy
        print 'p-value:', stats.chisqprob(CS,2)

def polarFind(theta, nbf, alpha):
    global MLE, funEval
    nbf.setParm(MLE)
    r0 = 1
    pe = funEval(r0,theta,nbf,alpha)
    print 'pe', pe
    while pe > 0.0:
        r0 = r0*5
        pe = funEval(r0,theta,nbf,alpha)
        print 'peI', pe
        assert(r0 < 500000)
    print 'f(a)', funEval(0,theta,nbf,alpha)
    print 'f(b)', funEval(r0,theta,nbf,alpha)
    r = optimize.brentq(funEval,0,r0,(theta,nbf,alpha))
    return r

def vectorFind(theta,nbf,alpha):
    global MLE
    r = polarFind(theta,nbf,alpha)
    v = h.Vector(2)
    v.x[0] = MLE[0] + r*math.cos(theta)
    v.x[1] = MLE[1] + r*math.sin(theta)
    return v

def polarInitial(r0,theta,nbf,alpha):
    global MLE, funEval
    nbf.setParm(MLE)
    z = funEval(r0,theta,nbf,alpha)
    if z > 0.0:
        a = r0
        b = r0*1.1
        while funEval(b,theta,nbf,alpha) > 0:
            b= b*1.1
            assert(b<1e9)
    elif z<0.0:
        b = r0
        a = r0*0.9
        while funEval(a,theta,nbf,alpha) < 0:
            a = a*0.9
            assert(a>1e-9)
    r = optimize.brentq(funEval,a,b,(theta,nbf,alpha))
    return r


def polarContinue(theta0, nbf, alpha, points):
    r = polarFind(theta0, nbf, alpha)
    rt = [[r,theta0]]
    for i in range(points):
        theta = theta0 + 2*math.pi*(i+1)/points
        return2MLE(nbf)
        print 'point #', i
        r = polarInitial(r,theta,nbf,alpha)
        # r = polarFind(theta, nbf, alpha)
        rt.append([r,theta])
    return rt

def polar2xy(rt,origin=None):
    global MLE
    xy = []
    if origin == None:
        M = MLE
    else:
        M = origin
        xy.append([M[0],M[1]])
    for rtheta in rt:
        xy.append([M[0] + rtheta[0]*math.cos(rtheta[1]), M[1] + rtheta[0]*math.sin(rtheta[1])])
    return xy

def log2linear(List):
    linxy = []
    for xy in List:
        linxy.append([math.exp(xy[0]), math.exp(xy[1])])
    return linxy

def testxy(nbf, List):
    for xy in List:
        p = h.Vector(2)
        p.x[0] = xy[0]
        p.x[1] = xy[1]
        nbf.setParm(p)
        pv = pValue(nbf,p)
        print 'pv', pv

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

funEval = hessEval
