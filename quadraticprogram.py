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

    def setGUIConstraints(self,geq0,leq1,sumto1):
        dim = len(geq0)
        nsums = len(sumto1)
        assert(dim == len(leq1))
        for j in range(nsums):
            assert(dim == len(sumto1[j]))
        geq0list = []
        leq1list = []
        sumto1list = []
        for j in range(nsums):
            sumto1list.append([])
        for i in range(dim):
            if geq0[i].x:
                geq0list.append(i)
            if leq1[i].x:
                leq1list.append(i)
            for j in range(nsums):
                if sumto1[j][i].x:
                    sumto1list[j].append(i)
        self.geq = (len(geq0list) > 0)
        self.leq = (len(leq1list) > 0)
        self.ineq = self.leq | self.geq
        self.eq = False
        for j in range(nsums):
            self.eq = self.eq | (len(sumto1list[j]) > 1)
        self.anyConstraints = self.eq | self.ineq
        if not self.eq:
            sumto1list = None
        self.setNeuronConstraints(dim,geq0list,leq1list,sumto1list)
        # print 'geq0list', geq0list
        # print 'leq1list', leq1list
        # print 'sumto1list', sumto1list
        # print 'Eq', self.eq
        # print 'Ineq', self.ineq
        # print 'Any', self.anyConstraints

    def setNeuronConstraints(self,dim,geq0,leq1,sum1):
        zeroRow = numpy.matlib.zeros((1,dim))
        start = True
        i = 0
        if sum1 == None:
            Deq = None
            eqd = None
        else:
            for sumlist in sum1:
                if len(sumlist) > 1:  # Must have at least to elements
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
        if not self.geq:
            Dgeq = None
            geqd = None
        else:
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
        if not self.leq:
            Dleq = None
            leqd = None
        else:
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

    def setConstraints(self,Deq,eqd,Dleq,leqd,Dgeq,geqd):
        if self.eq:
            self.A = cvxopt.matrix(Deq)
            self.b = cvxopt.matrix(eqd)
            # print 'A', self.A
            # print 'b', self.b
        if not self.leq and not self.leq:
            G = None
            h = None
        elif not self.leq:
            G = -Dgeq
            h = -geqd
        elif not self.geq:
            G = Dleq
            h = leqd
        else:
            G = numpy.concatenate((Dleq,-Dgeq))
            h = numpy.concatenate((leqd,-geqd))
        if not G == None:
            self.G = cvxopt.matrix(G)
            self.h = cvxopt.matrix(h)

    def setObjective(self,P,q):
        self.P = cvxopt.matrix(P)
        self.q = cvxopt.matrix(q)

    def project(self,m,PI=None):
        if PI == None:
            dim = numpy.matrix(self.G).shape[1]
            PI = numpy.matlib.eye(dim)
        self.setObjective(PI, -PI*m)
        return self.solve()

    def solve(self):
        if self.eq:
            out = solvers.qp(self.P,self.q,self.G,self.h,self.A,self.b)
        else:
            out = solvers.qp(self.P,self.q,self.G,self.h)
        return numpy.matrix(out['x'])