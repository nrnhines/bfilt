import cvxopt
import cvxopt.solvers as solvers
import numpy
import numpy.matlib

class QuadraticProgram(object):
    def __init__(self):
        self.P = None
        self.q = None
        self.G = None
        self.h = None
        self.A = None
        self.b = None
        self.eq = False
        solvers.options['show_progress'] = 0

    def setConstraints(self,Deq,eqd,Dleq,leqd,Dgeq,geqd):
        if Deq == None or eqd == None:
            self.eq = False
        else:
            self.eq = True
            self.A = cvxopt.matrix(Deq)
            self.b = cvxopt.matrix(eqd)
        G = numpy.concatenate((Dleq,-Dgeq))
        h = numpy.concatenate((leqd,-geqd))
        self.G = cvxopt.matrix(G)
        self.h = cvxopt.matrix(h)

    def setNeuronConstraints(self,dim,geq0,leq1,sum1):
        zeroRow = numpy.matlib.zeros((1,dim))
        start = True
        i = 0
        if sum1 == None:
            Deq = None
            eqd = None
        else:
            for sumlist in sum1:
                if start:
                    Deq = zeroRow.copy()
                    eqd = numpy.matrix([[1.0]])
                    start = False
                else:
                    Deq = numpy.concatenate((Deq,zeroRow.copy()))
                    eqd = numpy.concatenate((eqd,numpy.matrix([[1.0]])))
                for sumindex in sumlist:
                    Deq[i,sumindex] = 1.0
                i +=1
            start = True
            i = 0
        for component in geq0:
            if start:
                Dgeq = zeroRow.copy()
                geqd = numpy.matrix([[0]])
                start = False
            else:
                Dgeq = numpy.concatenate((Dgeq,zeroRow.copy()))
                geqd = numpy.concatenate((geqd,numpy.matrix([[0.0]])))
            Dgeq[i,component] = 1.0
            i += 1
        start = True
        i = 0
        for component in leq1:
            if start:
                Dleq = zeroRow.copy()
                leqd = numpy.matrix([[1]])
                start = False
            else:
                Dleq = numpy.concatenate((Dleq, zeroRow.copy()))
                leqd = numpy.concatenate((leqd,numpy.matrix([[1.0]])))
            Dleq[i,component] = 1.0
            i += 1
        self.setConstraints(Deq,eqd,Dleq,leqd,Dgeq,geqd)

    def setObjective(self,P,q):
        self.P = cvxopt.matrix(P)
        self.q = cvxopt.matrix(q)
    def solve(self):
        if self.eq:
            out = solvers.qp(self.P,self.q,self.G,self.h,self.A,self.b)
        else:
            out = solvers.qp(self.P,self.q,self.G,self.h)
        return numpy.matrix(out['x'])