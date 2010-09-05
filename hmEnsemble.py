import copy
import math
import numpy
import numpy.matlib

class Ensemble(object):
    def __init__(self,nchannels,pconfig,output_config,smallQ):
        self.nchannels = nchannels
        self.pconfig = pconfig
        self.output_config = output_config
        self.smallQ = smallQ
        self.nconfig = len(self.pconfig)
        self.enum = self.enumerate(self.nchannels, self.nconfig)
        self.nstates = len(self.enum)
        self.index = self.reverseEnumerate(self.nchannels, self.nconfig)
        self.pstates = self.makeProb(self.enum, self.pconfig)
        self.output = self.makeOutput(self.enum,self.output_config)
        self.Q = self.makeQ(self.enum,self.index,self.smallQ)

    def recursiveLookup(self,Ensem,n):
        if len(n) > 1:
            return self.recursiveLookup(Ensem[n[0]],n[1:])
        else:
            return Ensem[n[0]]

    def lookup(self,Ensem,n):
        return self.recursiveLookup(Ensem, n[:-1])  # last element in sequence is redundant and must be removed

    def makeQ(self,enum,index,smallQ):
        Q = numpy.matlib.zeros([len(enum), len(enum)])
        for e in enum:
            eIndex = self.lookup(index,e)
            for i in range(len(e)):
                if e[i] > 0:
                    for j in range(len(e)):
                        if not j == i:
                            f = copy.deepcopy(e)
                            f[i] -= 1   # New has one less in ith config
                            f[j] += 1  # New has one more in jth config, Note: j instead of i
                            fIndex= self.lookup(index,f)
                            newRate = smallQ[i,j]*e[i]
                            Q[eIndex,fIndex] = newRate
        # NOW DO i==j
        for i in range(len(enum)):
            sum = 0.0
            for j in range(len(enum)):
                sum += Q[i,j]
            Q[i,i] = -sum
        return Q

    def enumerate(self, nchannels, nconfig):  #Enumerate the possible states of an emsemble
        if nconfig > 1:
            enum = []
            for n0 in range(nchannels+1):
                enum0 = self.enumerate(nchannels-n0,nconfig-1) #recursive defn
                for e0 in enum0:
                    e0.append(n0)
                    enum.append(e0)
            return enum
        else:
            return [[nchannels]]

    def emptyEnsemble(self, nchannels, nconfig):   #Return an empty ensemble data structure
        if nconfig > 1:
            emp = []
            for n0 in range(nchannels+1):
                emp0 = self.emptyEnsemble(nchannels-n0,nconfig-1) # recursive defn
                emp.append(emp0)
            return emp
        else:
            return None

    def reverseEnumerate(self, nchannels, nstates):
        ind = self.emptyEnsemble(nchannels, nstates) # Data structure for reverse enumeration
        # e.g. Index[a,b,c] has number in enumeration corresponding to sequence [a,b,c]
        for i in range(self.nstates):
            self.assign(ind, self.enum[i], i)
        return ind

    def makeProb(self, enum, pconfig):
        nconfig = len(pconfig)
        nstates = len(enum)
        pstates = []
        for i in range(self.nstates):
            pstates.append(self.multiTerm(enum[i],pconfig))
        return pstates
        # pstatesIndex = self.emptyEnsemble(self.nchannels, nconfig)
        # for i in range(self.nstates):
            # self.assign(pstatesIndex, enum[i], self.multiTerm(enum[i],pconfig))
        # return (pstates, pstatesIndex)

    def recursiveAssign(self, Ensem, n, i):
        if len(n) > 1:
            self.recursiveAssign(Ensem[n[0]],n[1:],i)  # We have to make our way through the Ensem data structure
        else:
            Ensem[n[0]] = i

    def assign(self, Ensem, n, i):
        self.recursiveAssign(Ensem, n[:-1], i)  # last element in sequence is redundant and must be removed

    def multiTerm(self,n,p):
        # Used to derive ensemble probabilities from single probabilties
        # n is the list of numbers of channels in each state
        # p is a list of initial probabilities for each single channel state
        self.multiCheck(n,p)
        term = math.log(math.factorial(sum(n)))
        for i in range(len(n)):
            term -= math.log(math.factorial(n[i]))
        for i in range(len(p)):
            if p[i] > 0.0:
                term += n[i]*math.log(p[i])
            elif n[i] > 0:  # log is -infinity hence return is 0, unless n[i] is also zero, then log is 0 no contribution to answer
                return 0.0
        return math.exp(term)

    def multiCheck(self,n,p):
        # See comments for multiTerm(n,p) for defn of n and p, above.
        tol = 1e-7
        assert(len(n) == len(p))
        assert(min(p) >= 0.0)
        assert(min(n) >= 0)
        assert(math.fabs(sum(p) - 1.0) < tol)

    def makeOutput(self,enum,out):
        output = []
        for e in enum:
            o = 0.0
            for j in range(len(out)):
                o += out[j]*e[j]
            output.append(o)
        return output

