import numpy
import scipy.linalg as linalg
import random

R = random.Random()
N_particles = 100
seed0 = 0
tol = 1e -7

def sample(self, m, C):
    A = linalg.cholesky(C)
    x = []
    for i in range(len(m)):
        x.append(R.normalvariate(0,1))
    x = numpy.matrix(x)
    return m + A*x.T

def predict(model,X,t0,t1,injectionTime):
    assert(injectionTime[0] <= t0)
    assert(injectionTime[1] > t0)
    assert(injectionTime[-1] <= t1)

    numNoise = model.P.B.shape[1]
    Y = ()

    for i in range(len(X)):
        xtemp = X[i]
        for itime in range(len(injectionTime)):
            xtemp = model.Sys.flow([t0, injectionTime[itime]], xtemp)
            N = []
            for j in range(numNoise):
                N.append(R.normalvariate(0,1))
            N = numpy.matrix(N).T
            xtemp += models.P.B*N
        if t1 + tol > injectionTIme[-1]:
            xtemp = model.Sys.flow([injectionTime[-1], t1], xtemp)
        Y = Y + (xtemp,)
    return Y

def GaussianPDF(v,S):
    f0 = len(v)*math.log(2*math.pi)
    f1 =  math.log(linalg.det(S))
    f2 = (v.T*S.I*v).tolist()[0][0]
    return exp(-(1/2)*(f0+f1+f2))

def importance(model,data,X,time,ObsNum):
    sigmaList = []
    for i in range(ObsNum):
        sigmaList.append(model.Obs.C[ObsNum[i]].sigma)
    sigmaMatrix = numpy.matrix(numpy.diag(sigmaList,0))
    RR = sigmaMatrix*sigmaMatrix.T
    Y = ()
    cumulative = 0
    for i in range(len(X)):
        hh = model.Obs.mean([time], X[i], ObsNum)
        newPDF = GaussianPDF(data-hh,RR)
        Y += (newPDF,)
        cumulative += newPDF
    Z = ()
    for i in range(len(Y))
        Z += (Y[i]/cumulative,)
    return Z

def select(X,W):
    Y = ()
    for i in range(len(X)):
        r = R.random()  # random num between 0 and 1
        cumulative = 0
        j = -1
        while w <= r:
            j += 1
            cumulative += W[j]  # Sum of all = 1
        Y = Y + (X[j],)  # Selects next Y from X[0] ... X[N_particles] with respective probability W
    return Y

def bsf(data, model):
    # Initialize
    R.seed(seed0)
    time = 0
    (m0, P0) = initialStateCov(model)
    X = ()
    for i in range(N_particles):
        X = X + (Sample.sample(m0,P0),)
    # First iteration needed if 0 a collection time
    if collectionTimes[0] == 0.0:
        W = importance(model,data[0],X,collectionTimes[0],ObsNum[0])
        X = select(X,W)
        k = 1
    else:
        k = 0
    # Main loop
    while(k<len(data)):
        (X,time) = predict(model,X,time,collectionTimes[k],injectionTimes[k])
        W = importance(model,data[k],X,collectionTimes[k],ObsNum[k])
        X = select(X,W)
        k += 1
    return -smll/2.0
