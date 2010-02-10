import numpy
import scipy.linalg as linalg
import random

class Sample(object):
    R = random.Random()

    def seed(self, i):
        self.R.seed(i)

    def sample(self, m, C):
        A = linalg.cholesky(C)
        x = [0]*len(m)
        for i in range(len(m)):
            x[i] = self.R.normalvariate(0,1)
        x = numpy.matrix(x)
        return m + A*x.T