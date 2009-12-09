import math
from myscipy import linalg
import numpy
import fitglobals

def initialErrorBars():
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
    
def measCovJacs(time,ObsNum,m,P):
# Returns the Measurement and Covariance plus Jacobians
    hh = model.Obs.mean(time, m, ObsNum)
    H = model.Obs.Dstate(time, m, ObsNum)
    V = model.Obs.Dnoise(time, m, ObsNum)
    S = H*P*H.T + V*V.T
    return (hh,S,H,V)

def update(data,time,ObsNum,mb,Pb)
    (hh,S,H,V) = measCovJacs(time,ObsNum,m,P):
    e = data - hh
    K = Pb*H.T*S.I
    P = Pb - K*S*K.T
    m = mb + K*e

def ekf(data, model):
    # Initialize
    global Ecenter, Ewidth, Etime
    (Etime,Ecenter,Ewidth) = initialErrorBars()
    (m0, P0) = initialStateCov(model)
    k = 0
    smll = 0

    # main loop
    while(k<len(data)):
        if collectionTime[0] == 0:           