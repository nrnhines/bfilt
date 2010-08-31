import math
from myscipy import linalg
import numpy
import fitglobals
import HHBounds
import svd
import quadraticprogram
import copy
from neuron import h

def initializeErrorBars(Obs,Sys):
    global saveErrorBars, Etime, Ecenter, Ewidth
    saveErrorBars = True
    Etime = []
    Ecenter = []
    Ewidth = []
    for i in range(Obs.D):
        Ecenter.append([])
        Ewidth.append([])
    return (Etime,Ecenter,Ewidth)

def modelMeasurement(Obs,time,ObsNum,m,P):
    # Returns the measurement and covariance plus Jacobians
    hh = Obs.mean([time], m, ObsNum)
    H = Obs.Dstate([time], m, ObsNum)
    V = Obs.Dnoise([time], m, ObsNum)
    S = H*P*H.T + V*V.T
    return (hh,S,H,V)

def saveData(Obs,time,m,P):
    global saveErrorBars
    if saveErrorBars:
        ObsNum = range(Obs.D)
        hh = Obs.mean([time], m, ObsNum)
        H = Obs.Dstate([time], m, ObsNum)
        V = Obs.Dnoise([time], m, ObsNum)
        S = H*P*H.T + V*V.T
        global Etime, Ecenter, Ewidth
        Etime.append(time)
        for iObs in range(Obs.D):
            Ecenter[iObs].append(hh[iObs,0])
            Ewidth[iObs].append(math.sqrt(S[iObs,iObs]))

def dumpFunnels():
    global Etime, Ecenter, Ewidth
    f = open('funnels.txt','w')
    for i in range(len(Etime)):
        f.write('{0:f} {1:f} {2:f}\n'.format(Etime[i],Ecenter[0][i],Ewidth[0][i]))
    f.close()

def update(Obs,data,time,ObsNum,mb,Pb,bounds):
    (hh,S,H,V) = modelMeasurement(Obs,time,ObsNum,mb,Pb)
    saveData(Obs,time,mb,Pb)
    e = data - hh
    K = Pb*H.T*S.I
    P = Pb - K*S*K.T
    m = mb + K*e
    if QP.anyConstraints:
        PI = P.I
        QP.setObjective(PI,-PI*m)
        m = QP.solve()
    saveData(Obs,time,m,P)  # Saves error bars
    return (m,P,e,S)

def predict(Eve,Sys,m,P,t0,t1,injectionTimes):
    (mbs,times,oneStepDsF) = oneStepDFlowTable(t0,t1,injectionTimes,m,Sys)  # State times and Jacobians for one step intervals
    Bs = injectionEffectsList(injectionTimes, Eve)
    Pbr = PbTable(P,oneStepDsF,Bs)
    assert(len(Pbr) == len(mbs))
    assert(len(Pbr) == len(times))
    for i in range(len(Pbr)):  # starts at one because we have already handled initial point
        saveData(Eve.Obs,times[i],mbs[i],Pbr[i])
    return (mbs[-1],Pbr[-1],t1)

def minusTwiceLogGaussianPDF(v,S):
    f0 = len(v)*math.log(2*math.pi)
    f1 =  math.log(linalg.det(S))
    f2 = (v.T*S.I*v).tolist()[0][0]
    return (f0+f1+f2)

def ekf(data, Eve, Sys, DLikeDt_hvec = None):
    # Initialize
    if DLikeDt_hvec != None:
        DLikeDt_hvec.resize(0)
    smll = 0.0
    time = 0.0
    initializeErrorBars(Eve.Obs, Sys)
    collectionTimes = Eve.collectionTimes
    injectionTimes = Eve.injectionTimes
    print 'InjectionTable', injectionTimes
    bounds = HHBounds.bounds # model.stateBoundaries
    ObsNum = Eve.ObsNum
    (m0, P0) = initialStateCov(Eve.Sto,Sys)

    # Main loop
    if collectionTimes[0] == 0.0:
        (m,P,e,S) = update(Eve.Obs,data[0],collectionTimes[0],ObsNum[0],m0,P0,bounds)
        mll = minusTwiceLogGaussianPDF(e,S)
        if DLikeDt_hvec != None:
            DLikeDt_hvec.append(mll*0.5)
        smll += mll
        k = 1
    else:
        (m,P) = (m0,P0)
        k = 0
    while(k<len(data)):
        (mb,Pb,time) = predict(Eve,Sys,m,P,time,collectionTimes[k],injectionTimes[k])
        (m,P,e,S) = update(Eve.Obs,data[k],collectionTimes[k],ObsNum[k],mb,Pb,bounds)
        mll = minusTwiceLogGaussianPDF(e,S)
        if DLikeDt_hvec != None:
            DLikeDt_hvec.append(mll*0.5)
        smll += mll
        k += 1
    return -smll/2.0
