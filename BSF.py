import numpy
import scipy.linalg as linalg
import random
import math

R = random.Random()
N_particles = 100
seed0 = 0
tol = 1e-7

def initialStateCov(model):
    m0 = model.Initial
    P0 = model.P.InitialCov
    return (m0, P0)

def sample(m, C):
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
            xtemp += model.P.B*N
        if t1 + tol > injectionTime[-1]:
            xtemp = model.Sys.flow([injectionTime[-1], t1], xtemp)
        Y = Y + (xtemp,)
    return (Y, t1)

def GaussianPDF(v,S):
    f0 = len(v)*math.log(2*math.pi)
    f1 =  math.log(linalg.det(S))
    f2 = (v.T*S.I*v).tolist()[0][0]
    return math.exp(-(1/2)*(f0+f1+f2))

def initializeCloud(model):
    global Ctime, Cpoint
    Ctime = []
    Cpoint = []
    for i in range(model.Obs.D):
        Ctime.append([])
        Cpoint.append([])
    return (Ctime,Cpoint)

def saveCloud(time,hh,ObsNum):
    global Ctime, Cpoint
    for iObs in range(len(ObsNum)):
        Cpoint[ObsNum[iObs]].append(hh[iObs,0])
        Ctime[ObsNum[iObs]].append(time)

def initializeCloud2(model):
    global C2time, C2point
    C2time = []
    C2point = []
    for i in range(model.Obs.D):
        C2time.append([])
        C2point.append([])
    return (C2time,C2point)

def saveCloud2(time,H,ObsNum):
    global C2time, C2point
    for hh in H:
        for iObs in range(len(ObsNum)):
            C2point[ObsNum[iObs]].append(hh[iObs,0])
            C2time[ObsNum[iObs]].append(time)

def importance(model,data,X,time,ObsNum):
    sigmaList = []
    for i in range(len(ObsNum)):
        sigmaList.append(model.Obs.C[ObsNum[i]].sigma)
    sigmaMatrix = numpy.matrix(numpy.diag(sigmaList,0))
    RR = sigmaMatrix*sigmaMatrix.T
    Y = ()
    H = ()
    cumulative = 0
    for i in range(len(X)):
        hh = model.Obs.mean([time], X[i], ObsNum)
        newPDF = GaussianPDF(data-hh,RR)
        Y += (newPDF,)
        H += (hh,)  # H is measurement under X
        cumulative += newPDF
        saveCloud(time,hh,ObsNum)  # Could also include X[i] to save state
    Z = ()
    for i in range(len(Y)):
        Z += (Y[i]/cumulative,)
    return (Z,H)

def select(X,H,W):  # H is meaurement under X
    Y = ()
    G = ()
    for i in range(len(X)):
        r = R.random()  # random num between 0 and 1
        cumulative = 0
        j = -1
        while cumulative <= r:
            j += 1
            cumulative += W[j]  # Sum of all equals 1
        Y = Y + (X[j],)  # Selects next Y from X[0] ... X[N_particles] with respective probability W
        G = G + (H[j],)  # Select corresponding H
    return (Y,G)

def bsf(data, model):
    # Initialize
    initializeCloud(model)
    initializeCloud2(model)
    R.seed(seed0)
    time = 0
    collectionTimes = model.collectionTimes
    injectionTimes = model.injectionTimes
    ObsNum = model.ObsNum
    (m0, P0) = initialStateCov(model)
    X = ()
    saveX = []  # If saving state
    for i in range(N_particles):
        X = X + (sample(m0,P0),)
    # First iteration needed if 0 a collection time
    if collectionTimes[0] == 0.0:
        (W,H) = importance(model,data[0],X,collectionTimes[0],ObsNum[0])
        # saveX.append(X)    # If saving state after
        (X,H) = select(X,H,W)
        saveCloud2(collectionTimes[0],H,ObsNum[0])
        k = 1
    else:
        k = 0
    # Main loop
    while(k<len(data)):
        (X,time) = predict(model,X,time,collectionTimes[k],injectionTimes[k])
        (W,H) = importance(model,data[k],X,collectionTimes[k],ObsNum[k])
        (X,H) = select(X,H,W)
        saveCloud2(collectionTimes[k],H,ObsNum[k])
        k += 1
    return
