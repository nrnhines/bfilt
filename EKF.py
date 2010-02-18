import math
from myscipy import linalg
import numpy
import fitglobals
import HHBounds

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

def updateInBounds(K,e,mb,bounds):
# Boundaries are B*x >= b where B=bounds[i][0], B is a row-matrix, b=bounds[i][1], x=State
# The default update (m0 = mb + K*e) is adjusted to m0 = mb + alpha*K*e
    alpha = 1  # default update assuming its in bounds
    tolfactor = 1  # tolfactor reduces alpha more to prevent numerical issues
    tolbound = 1e-7 # tolerance for evaluating boudnaries
    Ke = K*e;  # need this more than once
    if Ke.T*Ke != 0: # Trap case where no update is needed, would break code
        for i in range(len(bounds)):
            # print 'i', i
            # print 'bounds[i][1]', bounds[i][1]
            # print 'bounds[i][0]', bounds[i][0]
            # print 'mb', mb
            # print 'Ke', Ke
            assert(bounds[i][0]*mb + tolbound >= bounds[i][1])
            # solves for alpha such that update on boundary
            newalpha = (bounds[i][1]-bounds[i][0]*mb)/(bounds[i][0]*Ke)
            #print 'newalpha', newalpha
            # reassigns alpha if a smaller update is required for this boundary
            if alpha > newalpha and newalpha > 0:
                alpha = newalpha
                tolfactor = 0.99999 # reduce final alpha by this amount
        assert(alpha >= 0)
    # print 'alpha', alpha
    # print 'tolfactor', tolfactor
    # print 'Ke', Ke
    return Ke*alpha*tolfactor

def update(model,data,time,ObsNum,mb,Pb,bounds):
    (hh,S,H,V) = modelMeasurement(model,time,ObsNum,mb,Pb)
    e = data - hh
    K = Pb*H.T*S.I
    P = Pb - K*S*K.T
    m = mb + updateInBounds(K,e,mb,bounds)
    modelMeasurement(model,time,ObsNum,m,P)  # Saves error bars
    return (m,P,e,S)

def predict(model,m,P,t0,t1,injectionTime):
    tol = 1e-7
    assert(injectionTime[0] <= t0 + tol)
    assert(injectionTime[1] > t0)
    assert(injectionTime[-1] <= t1 + tol)
    mb = m
    tStart = t0
    identityMatrixSizeOfState = numpy.eye(len(m))
    As = []
    Bs = []
    mbs = []
    for i in range(1,len(injectionTime)):
        (mb, A, B, tStart) = model.Sys.flowJac(tStart, [injectionTime[i-1],injectionTime[i]],mb)
        As.append(A)  # A's are Jacobians of above flows
        Bs.append(B)  # Typically all B's same matrix scaled by sqrt(dt)
        mbs.append(mb)
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
        Am_temp = Am*As[0]
        Pb_temp = Wm*Wm.T + Am_temp*P*Am_temp.T
        modelMeasurement(model,injectionTime[i+1],[0],mbs[i],Pb_temp)  # Saves error bars ONLY ObsNum = 0
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
    bounds = HHBounds.bounds # model.stateBoundaries
    ObsNum = model.ObsNum
    (m0, P0) = initialStateCov(model)

    # Main loop
    if collectionTimes[0] == 0.0:
        (m,P,e,S) = update(model,data[0],collectionTimes[0],ObsNum[0],m0,P0,bounds)
        k = 1
    else:
        (m,P) = (m0,P0)
        k = 0
    while(k<len(data)):
        (mb,Pb,time) = predict(model,m,P,time,collectionTimes[k],injectionTimes[k])
        (m,P,e,S) = update(model,data[k],collectionTimes[k],ObsNum[k],mb,Pb,bounds)
        smll += minusTwiceLogGaussianPDF(e,S)
        k += 1
    return -smll/2.0

