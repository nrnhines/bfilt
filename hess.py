import numpy
import numpy.matlib
import numdifftools as nd
import math
from neuron import h

def testfun(x,y):
    return 10*((x+y)**2) + 30*((x-y)**2)

def evalfun(p, N):
    saveParm = N.getParm()
    newParm = N.getParm()
    for i in range(len(p)):
        newParm.x[i] = p[i]
    N.setParm(newParm)
    L = N.likelihood()
    N.setParm(saveParm)
    return L

def accurateHessian(parm, N):
    parmList = []
    for i in range(int(parm.size())):
        parmList.append(parm[i])
    # Create (likelihood) inline function
    LamFun = lambda p: evalfun(p, N)
    # Create Hessian (of likelihood) inline function
    HessFun = nd.Hessian(LamFun)
    # Evaluate Hessian and return
    return numpy.matrix(HessFun(parmList))

def quickHessian(parm, N):
    print 'parm'
    parm.printf()
    sqrteps = math.sqrt(numpy.finfo(numpy.double).eps)
    nparms = int(parm.size())
    H = numpy.matlib.zeros([nparms,nparms])
    p = h.Vector()
    p.copy(parm)
    f = evalfun(p, N)
    p.printf()
    print 'f', f
    fi = []
    for i in range(nparms):
        p.x[i] += sqrteps
        fi.append(evalfun(p,N))
        p.printf()
        print 'fi', i, fi[i]
        p.x[i] -= sqrteps
    for i in range(nparms):
        for j in range(nparms):
            if i>j:
                H[i,j] = H[j,i]
            else:
                p.x[i] += sqrteps
                p.x[j] += sqrteps
                fij = evalfun(p,N)
                p.printf()
                print 'fij', i, j, fij
                p.x[i] -= sqrteps
                p.x[j] -= sqrteps
                H[i,j] = (fij - fi[i] - fi[j] + f)/(sqrteps**2.)
    return H