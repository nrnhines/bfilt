import math
from myscipy import linalg
import numpy
import fitglobals

def initializeErrorBars():
    global Etime, Ecenter, Ewidth
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
    
def modelMeasurement(time,ObsNum,m,P,saveErrorBars):
# Returns the measurement and movariance plus Jacobians
    hh = model.Obs.mean(time, m, ObsNum)
    H = model.Obs.Dstate(time, m, ObsNum)
    V = model.Obs.Dnoise(time, m, ObsNum)
    S = H*P*H.T + V*V.T
    
    if saveErrorBars:
        global Etime, Ecenter, Ewidth
        Etime.append(time0)
        for iObs in range(len(ObsNum)):
            Ecenter[ObsNum[iObs]].append(hh[iObs,0])
            Ewidth[ObsNum[iObs]].append(math.sqrt(S[iObs,iObs]))
    return (hh,S,H,V)

def update(data,time,ObsNum,mb,Pb,saveErrorBars)
    (hh,S,H,V) = modelMeasurement(time,ObsNum,mb,Pb,saveErrorBars):
    e = data - hh
    K = Pb*H.T*S.I
    P = Pb - K*S*K.T
    m = mb + K*e
    return (m,P,e,S)

def predict(m,P,t0,t1)
    # NEED TO WRITE THIS FUNCTION
    
def ekf(data, model):
    # Initialize
    initializeErrorBars()
    (m0, P0) = initialStateCov(model)
    k = 0
    smll = 0
    SaveErrorBars = True
    previousTime = 0
    
    # Main loop
    while(k<len(data)):
        if k == 0 and collectionTimes[0] == 0:
            (m,P,e,S) = update(data[0],collectionTimes[0],ObsNums[0],m0,P0,SaveErrorBars)
            modelMeasurement(collectionTimes[0],ObsNums[0],m,P,SaveErrorBars)
            k = 1
        else:
            (m,P) = (m0,P0)
        (mb,Pb) = predict(m,P,previousTime,collectionTimes[k])
        (m,P,e,S) = update(data[k],collectionTimes[k],ObsNums[k],mb,Pb,SaveErrorBars)
        modelMeasurement(collectionTimes[k],ObsNums[k],m,P,SaveErrorBars)
        previousTime = collectionTimes[k]
        k += 1
