import math
from myscipy import linalg
import numpy
import fitglobals

def initializeErrorBars(model):
    global saveErrorBars, Etime, Ecenter, Ewidth
    saveErrorBars = True
    Etime = []
    Ecenter = []
    Ewidth = []   
    for i in range(model.Obs.D):
        Ecenter.append([])
        Ewidth.append([])
    return (Etime,Ecenter,Ewidth)

def initialStateCov(model):
    m0 = model.Initial
    P0 = model.P.InitialCov
    return (m0, P0)
    
def modelMeasurement(model,time,ObsNum,m,P):
# Returns the measurement and covariance plus Jacobians
    hh = model.Obs.mean([time], m, ObsNum)
    H = model.Obs.Dstate([time], m, ObsNum)
    V = model.Obs.Dnoise([time], m, ObsNum)
    S = H*P*H.T + V*V.T
    
    global saveErrorBars
    if saveErrorBars:
        global Etime, Ecenter, Ewidth
        Etime.append(time)
        for iObs in range(len(ObsNum)):
            Ecenter[ObsNum[iObs]].append(hh[iObs,0])
            Ewidth[ObsNum[iObs]].append(math.sqrt(S[iObs,iObs]))
    return (hh,S,H,V)

def update(model,data,time,ObsNum,mb,Pb):
    (hh,S,H,V) = modelMeasurement(model,time,ObsNum,mb,Pb)
    e = data - hh
    K = Pb*H.T*S.I
    P = Pb - K*S*K.T
    m = mb + K*e
    modelMeasurement(model,time,ObsNum,m,P)  # Saves error bars
    return (m,P,e,S)

def predict(model,m,P,t0,t1,injectionTime):
    assert(injectionTime[0] <= t0)
    assert(injectionTime[1] > t0)
    assert(injectionTime[-1] <= t1)
    mb = m
    tStart = t0
    identityMatrixSizeOfState = numpy.eye(len(m))
    As = []
    Bs = []
    for i in range(1,len(injectionTime)):
        print 'itM', injectionTime[i-1]
        print 'itI', injectionTime[i]
        (mb, A, B, tStart) = model.Sys.flowJac(tStart, [injectionTime[i-1],injectionTime[i]],mb)
        As.append(A)  # A's are Jacobians of above flows
        Bs.append(B)  # Typically all B's same matrix scaled by sqrt(dt)
    if  t1 > injectionTime[-1]:
        (mb, A, B, tStart) = model.Sys.flowJac(tStart,[injectionTime[-1],t1],mb)
        As.append(A)
    else:
        As.append(identityMatrixSizeOfState)
    Am = identityMatrixSizeOfState
    for i in range(len(Bs)):
        Am = Am*As[-(i+1)]  # Composition of Jacobians is product
        New = Am*Bs[-(i+1)]
        if i == 0:
            Wm = New
        else:
            Wm = numpy.bmat('New Wm')
    Am = Am*As[0]
    Pb = Wm*Wm.T + Am*P*Am.T
    return (mb, Pb, t1)
    
def minusTwiceLogGaussianPDF(v,S):
    f0 = len(v)*math.log(2*math.pi)
    f1 =  math.log(linalg.det(S))
    f2 = (v.T*S.I*v).tolist()[0][0]
    return (f0+f1+f2)

def ekf(data, model):
    # Initialize
    smll = 0.0
    time = 0.0
    initializeErrorBars(model)
    collectionTimes = model.collectionTimes
    injectionTimes = model.injectionTimes
    ObsNum = model.ObsNum
    (m0, P0) = initialStateCov(model)
    
    # Main loop
    if collectionTimes[0] == 0.0:
        (m,P,e,S) = update(model,data[0],collectionTimes[0],ObsNum[0],m0,P0)
        k = 1
    else:
        (m,P) = (m0,P0)
        k = 0
    while(k<len(data)):
        (mb,Pb,time) = predict(model,m,P,time,collectionTimes[k],injectionTimes[k])
        (m,P,e,S) = update(model,data[k],collectionTimes[k],ObsNum[k],mb,Pb)
        smll += minusTwiceLogGaussianPDF(e,S)
        k += 1
    return -smll/2.0