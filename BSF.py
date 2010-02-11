import numpy
import scipy.linalg as linalg
import random

class gaussSample(object):
    R = random.Random()

    def seed(self, i):
        self.R.seed(i)

    def sample(self, m, C):
        A = linalg.cholesky(C)
        x = []
        for i in range(len(m)):
            x.append(self.R.normalvariate(0,1))
        x = numpy.matrix(x)
        return m + A*x.T


def bsf(data, model):
    # Constants
    N = 100
    seed0 = 0

    # Initialize
    Sample = gaussSample()
    Sample.seed(seed0)
    (m0, P0) = initialStateCov(model)
    X = ()
    for i in range(N):
        X = X + (Sample.sample(m0,P0),)

    #Main loop
    if collectionTimes[0] == 0.0:
        W = weights(model,data[0],X,collectionTimes[0],ObsNum[0])
        X = resample(X,W)
        k = 1
    else:
        k = 0
    while(k<len(data)):
        X = predict(model,X,time,collectionTimes[k],injectionTimes[k])
        W = weights(model,data[0],X,collectionTimes[0],ObsNum[0])
        X = resample(X,W)
        k += 1
    return -smll/2.0
