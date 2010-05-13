import math
from myscipy import linalg
import numpy
import fitglobals
import HHBounds
import svd

def constraintsOn(sum):
    global useConstraints
    if sum:
        (Deq,eqd) = equalityConstraints()
    else
        Deq = None
        eqd = None
    Dmatrix = numpy.matrix([[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0])
    d1matrix = numpy.matrix([[1.0],[1.0],[1.0]])
    d0matrix = numpy.matrix([[0.0],[0.0],[0.0]])
    Dgeq = Dmatrix.copy()
    geqd = d0matrix.copy()
    Dleq = Dmatrix.copy()
    leqd = d1matrix.copy()
    QP = QuadraticProgram()
    QP.setConstraints(Deq,eqd,Dleq,leqd,Dgeq,geqd)
    useConstraints = True


def constraintsOff():
    global useConstraints
    useConstraints = False

def initializeErrorBars(Obs,Sys):
    global saveErrorBars, Etime, Ecenter, Ewidth, Scenter, Swidth, Ps, ms
    saveErrorBars = True
    Etime = []
    Ecenter = []
    Ewidth = []
    Scenter = []
    Swidth = []
    Ps = []
    ms = []
    for i in range(Obs.D):
        Ecenter.append([])
        Ewidth.append([])
    for i in range(Sys.dim()):
        Scenter.append([])
        Swidth.append([])
    return (Etime,Ecenter,Ewidth)

def initialStateCov(Sto,Sys):
    m0 = Sys.Initial
    P0 = Sto.InitialCov
    return (m0, P0)

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
        global Etime, Ecenter, Ewidth, Scenter, Swidth, Ps, ms
        Etime.append(time)
        for iObs in range(Obs.D):
            Ecenter[iObs].append(hh[iObs,0])
            Ewidth[iObs].append(math.sqrt(S[iObs,iObs]))
        for s in range(P.shape[0]):
            Scenter[s].append(m[s])
            Swidth[s].append(math.sqrt(P[s,s]))
        Ps.append(P)
        ms.append(m)

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
            if (bounds[i][0]*mb + tolbound < bounds[i][1]):
                print 'i', i
                print 'bounds[i][1]', bounds[i][1]
                print 'bounds[i][0]', bounds[i][0]
                print 'mb', mb
                print 'tolbound', tolbound
                print 'bounds[i][0]*mb + tolbound', bounds[i][0]*mb + tolbound
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

def equalityConstraints():
    D = numpy.matrix([[0.0, 1.0, 1.0, 1.0]])
    d = numpy.matrix([[1]])
    return (D,d)

def addInequalitiesViolated(m,bounds,D,d):
    tol = 1e-7
    first = (D == None)
    allSatisfied = True
    for i in range(len(bounds)):
        if bounds[i][0]*m + tol < bounds[i][1]:
            allSatisfied = False
            Dnew = bounds[i][0]
            dnew = numpy.matrix(bounds[i][1])
            if first:
                D = Dnew
                d = dnew
                first = False
            else:
                # print 'd', d
                # print 'dnew', dnew
                D = numpy.bmat('D; Dnew')
                d = numpy.bmat('d; dnew')
    return (D,d,allSatisfied)

def project(m,P,D,d):
    # print 'project: m', m
    # print 'project: P', P
    # print 'project: D', D
    # print 'project: d', d
    mtilde = m - P*D.T*(D*P*D.T).I*(D*m - d)
    return mtilde

def update(Obs,data,time,ObsNum,mb,Pb,bounds):
    (hh,S,H,V) = modelMeasurement(Obs,time,ObsNum,mb,Pb)
    saveData(Obs,time,mb,Pb)
    e = data - hh
    K = Pb*H.T*S.I
    P = Pb - K*S*K.T
    # m = mb + updateInBounds(K,e,mb,bounds)
    m = mb + K*e
    global useConstraints
    if useConstraints:
        mold = m # REMOVE
        PI = P.I
        QP.setObjective(PI,-PI*m)
        m = QP.solve()
        # REMOVE TO END_IF
        (D,d) = equalityConstraints()
        (D,d,temp) = addInequalitiesViolated(mold,bounds,D,d)
        allSatisfied = (D==None)
        while not allSatisfied:
            mold = project(mold,P,D,d)
            (D,d,allSatisfied) = addInequalitiesViolated(mold, bounds, D, d)
        print "Constrained Diff", (mold - m).T*(mold - m)
    saveData(Obs,time,m,P)  # Saves error bars
    if fitglobals.debug:
        print 'New Pb'
        print 'Pb values', svd.svd(Pb.tolist())[1]
        print 'P values', svd.svd(P.tolist())[1]
    return (m,P,e,S)

def predict(Eve,Sys,m,P,t0,t1,injectionTime):
    tol = 1e-7
    assert(injectionTime[0] <= t0 + tol)
    assert(injectionTime[1] + tol > t0)
    assert(injectionTime[-1] <= t1 + tol)
    mb = m
    tStart = t0
    identityMatrixSizeOfState = numpy.eye(len(m))
    As = []
    Bs = []
    mbs = []
    for i in range(1,len(injectionTime)):
        (mb, A, tStart) = Sys.flowJac(tStart, [injectionTime[i-1],injectionTime[i]],mb)
        B = Eve.Sto.noiseJac([injectionTime[i-1],injectionTime[i]])
        As.append(A)  # A's are Jacobians of above flows
        Bs.append(B)  # Typically all B's same matrix scaled by sqrt(dt)
        mbs.append(mb)
    if  t1 > injectionTime[-1]:
        (mb, A, tStart) = Sys.flowJac(tStart,[injectionTime[-1],t1],mb)
        B = Eve.Sto.noiseJac([injectionTime[i-1],injectionTime[i]])
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
        saveData(Eve.Obs,injectionTime[i+1],mbs[i],Pb_temp)  # Saves error bars ONLY ObsNum = 0
    Am = Am*As[0]
    Pb = Wm*Wm.T + Am*P*Am.T
    return (mb, Pb, t1)

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
