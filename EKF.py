import math
from myscipy import linalg
import numpy
import fitglobals

def initializeErrorBars():
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
    
def modelMeasurement(time,ObsNum,m,P):
# Returns the measurement and movariance plus Jacobians
    hh = model.Obs.mean(time, m, ObsNum)
    H = model.Obs.Dstate(time, m, ObsNum)
    V = model.Obs.Dnoise(time, m, ObsNum)
    S = H*P*H.T + V*V.T
    
    global saveErrorBars
    if saveErrorBars:
        global Etime, Ecenter, Ewidth
        Etime.append(time0)
        for iObs in range(len(ObsNum)):
            Ecenter[ObsNum[iObs]].append(hh[iObs,0])
            Ewidth[ObsNum[iObs]].append(math.sqrt(S[iObs,iObs]))
    return (hh,S,H,V)

def update(data,time,ObsNum,mb,Pb):
    (hh,S,H,V) = modelMeasurement(time,ObsNum,mb,Pb):
    e = data - hh
    K = Pb*H.T*S.I
    P = Pb - K*S*K.T
    m = mb + K*e
    modelMeasurement(time,ObsNum,m,P)
    return (m,P,e,S)

def predict(m,P,t0,t1,injectionTimes):
    return (mb, Pb, t1)
    
def minusTwiceLogGaussianPDF(v,S):
    f0 = len(v)*math.log(2*math.pi)
    f1 =  math.log(linalg.det(S))
    f2 = (v.T*S.I*v).tolist()[0][0]
    return (f0+f1+f2)

def ekf(data, model):
    # Initialize
    smll = 0
    initializeErrorBars()
    (m0, P0) = initialStateCov(model)
    
    # Main loop
    if collectionTime[0] == 0:
        (m,P,e,S) = update(data[0],collectionTime[0],ObsNums[0],m0,P0)
        k = 1
    else:
        (m,P) = (m0,P0)
        k = 0
    while(k<len(data))        
        (mb,Pb,time) = predict(m,P,time,collectionTime[k],injectionTimes[k])
        (m,P,e,S) = update(data[k],collectionTime[k],ObsNum[k],mb,Pb,saveError)
        smll += minusTwiceLogGaussianPDF(e,S)
        k += 1
    return -smll/2.0