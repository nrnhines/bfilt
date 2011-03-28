from neuron import h
import testcr

def handle():
    return testcr.NrnBFiltHandle(testcr.MulRunFitHandle())

def quadfun(x,y,N=None):
    return (x+y)**2 + 0.5*(x-y)**2

def likeslice(x,y,tcr):
    saveParm = tcr.N.getParmVal()
    v = tcr.N.getParmVal()
    v.x[0] = x
    v.x[1] = y
    tcr.N.setParmVal(v)
    loglike = tcr.N.likelihood()
    tcr.N.setParmVal(saveParm)
    return loglike

def pslice(x,y,tcr):
    assert not (tcr == None)
    assert (tcr.lastrun >= 4)
    saveParm = tcr.N.getParmVal()
    v = tcr.N.getParmVal()
    v.x[0] = x
    v.x[1] = y
    tcr.N.setParmVal(v)
    minusPointLogLike = tcr.N.likelihood()
    tcr.N.setParmVal(saveParm)
    minusMaxLogLike = tcr.ml
    nParm = 2
    return tcr.get_pValue(minusPointLogLike,minusMaxLogLike,nParm)

def confslice(x,y,tcr):
    return 100*(1-pslice(x,y,tcr))

def evalpoint(x,y,nu=None, tcr=None, efun=likeslice):
    if tcr == None:
        assert False
        N = handle()
    else:
        N = tcr.N
    saveParm = N.getParmVal()
    if not nu == None:
        v = N.getParmVal()
        v.x[2] = nu
        N.setParmVal(v)
    z = efun(x,y,tcr)
    N.setParmVal(saveParm)
    return z

def eval2D(xx,yy,nu=None, tcr=None, efun=likeslice, onemain=False):
    if tcr == None:
        assert False
        N = handle()
    else:
        N = tcr.N
    saveParm = N.getParmVal()
    if not nu == None and not onemain:
        v = N.getParmVal()
        v.x[2] = nu
        N.setParmVal(v)
    zz = []
    for x in xx:
        zz.append([])
        for y in yy:
            zz[-1].append(efun(x,y,tcr))
            print 'x', x, 'y', y, 'z', zz[-1][-1]
    N.setParmVal(saveParm)
    return zz

def save2D(fname,xx,yy,nu=None,tcr=None,efun=confslice,onemain=False):
    zz = eval2D(xx,yy,nu,tcr,efun,onemain)
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