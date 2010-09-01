import copy
import math
import numpy
import numpy.matlib

def enumEnsemble(nchannels, nstates):  #Enumerate the possible states of an emsemble
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

def emptyEnsemble(nchannels, nstates):   #Return an empty ensemble data structure
    if nstates > 1:
        emp = []
        for n0 in range(nchannels+1):
            emp0 = emptyEnsemble(nchannels-n0,nstates-1) # recursive defn
            emp.append(emp0)
        return emp
    else:
        return None

def multiCheck(n,p):
    # See comments for multiTerm(n,p) below.
    tol = 1e-7
    assert(len(n) == len(p))
    assert(min(p) >= 0)
    assert(min(n) >= 0)
    assert(math.fabs(sum(p) - 1.0) < tol)

def multiTerm(n,p):
    # Used to derive ensemble probabilities from single probabilties
    # n is the list of numbers of channels in each state
    # p is a list of initial probabilities for each single channel state
    multiCheck(n,p)
    term = math.log(math.factorial(sum(n)))
    for i in range(len(n)):
        term -= math.log(math.factorial(n[i]))
    for i in range(len(p)):
        term += n[i]*math.log(p[i])
    return math.exp(term)

def recursiveAssignEnsemble(Ensem, n, i):
    if len(n) > 1:
        recursiveAssignEnsemble(Ensem[n[0]],n[1:],i)  # We have to make our way through the Ensem data structure
    else:
        Ensem[n[0]] = i

def assignEnsemble(Ensem, n, i):
    recursiveAssignEnsemble(Ensem, n[:-1], i)  # last element in sequence is redundant and must be removed

def recursiveLookup(Ensem,n):
    if len(n) > 1:
        return recursiveLookup(Ensem[n[0]],n[1:])
    else:
        return Ensem[n[0]]

def lookup(Ensem, n):
    return recursiveLookup(Ensem, n[:-1])  # last element in sequence is redundant and must be removed

def twoWayTable(nchannels, nstates):
    Enum = enumEnsemble(nchannels, nstates)  # Enumeration
    # Enum[i] is a list of channels in each state, e.g. [3,5,4]
    Index = emptyEnsemble(nchannels, nstates) # For reverse enumeration
    # e.g. Index[a,b,c] has number in enumeration corresponding to sequence [a,b,c]
    for i in range(len(Enum)):
        assignEnsemble(Index, Enum[i], i)
    return (Enum,Index)

def distribEnsemble(nchannels, nstates, prob):
    # From prob, the probability of a single channel being in each state
    # Derive eProb the ensemble probability
    (eEnum, eIndex) = twoWayTable(nchannels, nstates)
    eProb = emptyEnsemble(nchannels, nstates)
    for i in range(len(eEnum)):
        assignEnsemble(eProb, eEnum[i], multiTerm(eEnum[i],prob))
    return (eEnum, eIndex, eProb)

def assignRate(Rates,n,m,newRate):
    print 'hello'

def ratesEnsemble(eEnum,eIndex,rates,nchannels):
    eRates = numpy.matlib.zeros([nchannels, nchannels])
    for n in eEnum:
        for i in range(len(n)):
            if n[i] > 0:
                for j in range(len(n)):
                    if not j == i:
                        m = copy.deepcopy(n)
                        m[i] -= 1
                        m[j] += 1
                        newRate = rates[i,j]*n[i]
                        assignRate(eRates,n,m,newRate)
    #REMEMBER: COME BACK AND DO i==j

class Ensemble(object):
    def __init__(self,nchannels,nstates):
        (self.enum, self.index) = self.twoWayTable(nchannels, nstates)
        self.prob =
        self.Q
