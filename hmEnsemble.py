import copy
import math

def enumEnsemble(nchannels, nstates):
    if nstates > 1:
        enum = []
        for n0 in range(nchannels+1):
            enum0 = enumEnsemble(nchannels-n0,nstates-1) #recursive defn
            for e0 in enum0:
                e0.append(n0)
                enum.append(e0)
        return enum
    else:
        return [[nchannels]]

def emptyEnsemble(nchannels, nstates):
    if nstates > 1:
        emp = []
        for n0 in range(nchannels+1):
            emp0 = emptyEnsemble(nchannels-n0,nstates-1) # recursive defn
            emp.append(emp0)
        return emp
    else:
        return None

def multiCheck(n,p):
    tol = 1e-7
    assert(len(n) == len(p))
    assert(min(p) >= 0)
    assert(min(n) >= 0)
    assert(math.fabs(sum(p) - 1.0) < tol)

def multiTerm(n,p):
    multiCheck(n,p)
    term = math.log(sum(n))
    for i in range(len(n)):
        term -= math.log(n[i])
    for i in range(len(p)):
        term += n[i]*math.log(p[i])
    return math.exp(term)

def assignEnsemble(Ensem, n, i):
    if len(n) > 1:
        assignEnsemble(Ensem[n[0]],n[1:],i)
    else:
        Ensem[n[0]] = i

def twoWayTable(nchannels, nstates):
    Enum = enumEnsemble(nchannels, nstates)
    Index = emptyEnsemble(nchannels, nstates)
    for i in range(len(Enum)):
        assignEnsemble(Index, Enum[i][:-1], i)
    return (Enum,Index)

def distribEnsemble(nchannels, nstates, prob):
    eStates = enumEnsemble(nchannels, nstates)
    eProb = emptyEnsemble(nchannels, nstates)
