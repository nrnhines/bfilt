import cvxopt
import cvxopt.solvers as solvers
import numpy

class QuadraticProgram(object):
    def __init__(self):
        self.P = None
        self.q = None
        self.G = None
        self.h = None
        self.A = None
        self.b = None
        solvers.options['show_progress'] = 0
    def setConstraints(self,Deq,eqd,Dleq,leqd,Dgeq,geqd):
        self.A = cvxopt.matrix(Deq)
        self.b = cvxopt.matrix(eqd)
        G = numpy.concatenate((Dleq,-Dgeq))
        h = numpy.concatenate((leqd,-geqd))
        self.G = cvxopt.matrix(G)
        self.h = cvxopt.matrix(h)
    def setObjective(self,P,q):
        self.P = cvxopt.matrix(P)
        self.q = cvxopt.matrix(q)
    def solve(self):
        print 'self.b', self.b
        out = solvers.qp(self.P,self.q,self.G,self.h,self.A,self.b)
        return numpy.matrix(out['x'])