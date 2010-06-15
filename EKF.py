import math
from myscipy import linalg
import numpy
import fitglobals
import HHBounds
import svd
import quadraticprogram
import copy

def constraintsOn(geq0,leq1,sumto1):
    global QP
    QP = quadraticprogram.QuadraticProgram()
    QP.setGUIConstraints(geq0,leq1,sumto1)

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
    d = numpy.matrix([[1.0]])
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
                D = numpy.bmat('D; Dnew')
                d = numpy.bmat('d; dnew')
    return (D,d,allSatisfied)

def project(m,P,D,d):
    mtilde = m - P*D.T*(D*P*D.T).I*(D*m - d)
    return mtilde

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

    # OLD CONSTRAINTS CODE (goes before saveData above)
    #~ global useConstraints, QP, newConstraints, oldConstraints
    #~ if useConstraints:
        #~ mold = m
        #~ if newConstraints:
            #~ PI = P.I
            #~ QP.setObjective(PI,-PI*m)
            #~ m = QP.solve()
        #~ if oldConstraints:
            #~ (D,d) = equalityConstraints()
            #~ (D,d,temp) = addInequalitiesViolated(mold,bounds,D,d)
            #~ allSatisfied = (D==None)
            #~ while not allSatisfied:
                #~ mold = project(mold,P,D,d)
                #~ (D,d,allSatisfied) = addInequalitiesViolated(mold, bounds, D, d)
        #~ if newConstraints and oldConstraints:
            #~ print "Constrained Diff", (mold - m).T*(mold - m)
        #~ else:
            #~ m = mold

def oneStepDFlowTable(tStart, tFinal, injectionTimes, mb, Sys):
    tol = 1e-7
    assert(injectionTimes[0] <= tStart + tol) # injectionTime[0] is previous injection (possibly before tStart) (<= easier to satisfy with tol)
    assert(tStart + tol < injectionTimes[1])   # injectionTime[1] must be after tStart (first endpoint of step) (< harder to satisfy with tol)
    assert(injectionTimes[-1] <= tFinal + tol) # last injection time must be at or before tFinal (<= easier to satisfy with tol)
    oneStepDsF = []  # List of flow Jacobians to be returned by this function
    mbs = [mb]  # List of state vectors of flow (two end points to an interval so one more than in oneStepDsF) to be returned by this function
    times = [tStart]  # List of corresponding (to state) times to be returned by this function
    for i in range(1,len(injectionTimes)):
        mbs.append(mb)  # Before ith flow
        times.append(tStart) # Before ith flow
        (mb,A,tStart) = Sys.flowJac(tStart,[tStart,injectionTimes[i]],mb)  # overwrites mb and tStart, used next
        oneStepDsF.append(A) # Jacobian of ith flow (from point i to point i+1)
    mbs.append(mb)  # After flow to last injection time
    times.append(tStart) # After flow to last injection time
    if tFinal > (tStart + tol):  # final point if beyond last injection
        (mb,A,tFinal) = Sys.flowJac(tStart,[tStart,tFinal],mb)  # overwrites tFinal with same value, mb with different value
        mbs.append(mb) # After flow to final
        times.append(tFinal) # After flow to final
        oneStepDsF.append(A) # Jacobian of final flow
    return (mbs,times,oneStepDsF)

def multiStepDFlowTable(oneStepDsF):
    # Returns a table whose (i,j) element is the Jacobian of flow from i point i to point j as determined by oneStepDsF
    npoints = len(oneStepDsF) + 1
    nstates = oneStepDsF[0].shape[0]
    identityMatrixSizeOfState = numpy.matrix(numpy.eye(nstates))
    DsF = [] # Initialize table to be returned
    # Now fill (n) by (n) table with None.  Later DsF[i][j] will be Jacobian of flow wrt state from point i to point j
    for i in range(npoints):
        DsFi = [] # initialize ith row
        for j in range(npoints):
            DsFi.append(None)  # (i,j) element = None
        DsF.append(DsFi)  # add ith row
    # Now replace the None's with the elements that are needed
    # Leave as None DsF[i][j] with j<i (instead of what makes sense: DsF[j][i].I) because these values are not used
    for i in range(npoints):
        DsF[i][i] = identityMatrixSizeOfState  # Jacobians of trivial flows (time doesn't change) are identity matrix
    for i in range(npoints-1):
        DsF[i][i+1] = oneStepDsF[i]  # Fill in one-step Jacobians calculated earlier
    # Now compute the rest by composition
    n = npoints-1  # Only need computation for n=npoints-1 if you don't plot funnels, otherwise will repeat for other n
    for i in range(n-2,-1,-1):  # i.e. npoints-3,...,3,2,1,0, DsF[npoints-3][npoints-1] is first one of DsF[i][n] not defined
        DsF[i][n] = DsF[i+1][n]*DsF[i][i+1]
    # If plotting funnels repeat for the rest of the n's
    if True:  # Change this later so that its True only when plotting funnels
        for n in range(npoints-2,-1,-1): # i.e. npoints-2,npoints-1,...,3,2,1,0
            for i in range(n-2,-1,-1):
                DsF[i][n] = DsF[i+1][n]*DsF[i][i+1]
    return DsF

def injectionEffectsList(injectionTimes,Eve):
    Bs = []
    for i in range(len(injectionTimes)-1):
        Bs.append(Eve.Sto.noiseJac([injectionTimes[i],injectionTimes[i+1]]))
    return Bs

def lastTrivialStepAddedDFlowTable(Bs, DsF, injectionTimes, t1):
    # There may or may not have been an extra flow with no injection at the end
    # If there isn't, add a trivial flow to the end of DsF
    tol = 1e-7
    DsF2 = copy.deepcopy(DsF)
    nflows = len(DsF)  # DsF2 is a list of lists (of matricies)
    assert(nflows == len(DsF[0]))  # DsF2[i][j] should be square
    if (injectionTimes[-1] + tol < t1): # tol makes it harder to satisfy
        print 'Warning: No Injection at this Collection -- Code Untested'
    else:
        DsF2LastRow = []
        for i in range(nflows):
            DsF2[i].append(DsF[i][-1])    # DsF2[i][n+1] = DsF[i][n]
            DsF2LastRow.append(None)  # DsF2[n+1][i] = None
        DsF2LastRow.append(DsF[0][0]) # DsF2[n+1][n+1] = identity matrix size of state
        DsF2.append(DsF2LastRow)
    return DsF2

def DFlowWrtNoiseTable(Bs,DsF2):
    # DnF[i][j]: deriv wrt noise when injecting after ith injection interval then traversing the flow from the i+1 through jth flow point
    # i+1st flow point = ith injection point because it is the right endpoint of the ith injection interval, see picture:
    # Previous Inject ---<=--- Previous Collect ---FLOW0--->InjectionIntervalRightEndPoint0 ----FLOW1----> Inject1 ----FLOW2---> etc
    # There is always a last flow beyond the last injection, even if it is the trivial identity flow
    nrows = len(Bs)
    ncols = len(DsF2)
    assert(nrows + 2 == ncols)
    rowDnF = []  # append ncols None's to make a row of DnF
    for j in range(ncols):
        rowDnF.append(None)
    DnF = [] # append nrows rows of Nones to make DnF
    for i in range(nrows):
        DnF.append(rowDnF)
    for i in range(nrows):  # overwrite None's where possible
        for j in range(ncols):
            if DsF2[i+1][j] != None:
                DnF[i][j] = DsF2[i+1][j]*Bs[i]
    return DnF

def DFlowFromBeginWrtStateMatrixList(DsF2):
    # Ams[i] is Flow Jacobian when flowing from tStart to ith flow point
    Ams = []
    for i in range(0,len(DsF2)):
        Ams.append(DsF2[0][i])
    return Ams

def DFlowFromBeginWrtNoiseMatrixList(DnF):
    # Wms[i] is a matrix with several columns.
    # The jth column of Wms[i] is the derivative with respect to the jth injection interval after flowing from the j+1st to the ith flow point
    # For i=0, Wms[i] = None, because no injections before 0th flow point (right endpoint of first injection interval is flow point 1).
    Wms = [None]
    for i in range(1,len(DnF[0])):
        Wm = []
        WmEmpty = True
        for j in range(len(DnF)):
            if WmEmpty:
                if DnF[0][i] != None:
                    Wm = DnF[0][i]
                    WmEmpty = False
            else:
                if DnF[j][i] != None:
                    temp = DnF[j][i]
                    Wm = numpy.bmat('Wm temp')
        if WmEmpty:
            Wms.append(None)
        else:
            Wms.append(Wm)
    return Wms

def covarianceTable(Ams,Wms,P):
    Pbs = [P]  # initial covariance is P
    for i in range (1,len(Wms)): # starts at one because skipping initial point already handled
        Pbs.append(Wms[i]*Wms[i].T + Ams[i]*P*Ams[i].T)
    return Pbs

def predict(Eve,Sys,m,P,t0,t1,injectionTimes):
    (mbs,times,oneStepDsF) = oneStepDFlowTable(t0,t1,injectionTimes,m,Sys)  # State times and Jacobians for one step intervals
    DsF = multiStepDFlowTable(oneStepDsF)  # Table whose i,j element is Jacobian of flow from point i to point j
    Bs = injectionEffectsList(injectionTimes, Eve)
    DsF2 = lastTrivialStepAddedDFlowTable(Bs, DsF, injectionTimes, t1)
    DnF = DFlowWrtNoiseTable(Bs,DsF2)
    Ams = DFlowFromBeginWrtStateMatrixList(DsF2)
    Wms = DFlowFromBeginWrtNoiseMatrixList(DnF)
    Pbs = covarianceTable(Ams,Wms,P)
    for i in range(1,len(Pbs)):  # starts at one because we have already handled initial point
        saveData(Eve.Obs,times[i],mbs[i],Pbs[i])
    return (mbs[-1],Pbs[-1],t1)

    #~ identityMatrixSizeOfState = numpy.eye(len(m))
    #~ As = []
    #~ Bs = []
    #~ mbs = []
    #~ for i in range(1,len(injectionTime)):
        #~ (mb, A, tStart) = Sys.flowJac(tStart, injectionTime[i],mb)
        #~ B = Eve.Sto.noiseJac([injectionTime[i-1],injectionTime[i]])
        #~ # print 'Bs', B
        #~ As.append(A)  # A's are Jacobians of above flows
        #~ Bs.append(B)  # Typically all B's same matrix scaled by sqrt(dt)
        #~ mbs.append(mb)
    #~ if  t1 > injectionTime[-1]:
        #~ print 'WARNING: No Injection at this Collection'
        #~ (mb, A, tStart) = Sys.flowJac(tStart,[injectionTime[-1],t1],mb)
        #~ B = Eve.Sto.noiseJac([injectionTime[i-1],injectionTime[i]])
        #~ As.append(A)
    #~ else:
        #~ As.append(identityMatrixSizeOfState)
    #~ Am = identityMatrixSizeOfState
    #~ for i in range(len(Bs)):
        #~ Am = Am*As[-(i+1)]  # Composition of Jacobians is product
        #~ New = Am*Bs[-(i+1)]
        #~ if i == 0:
            #~ Wm = New
        #~ else:
            #~ Wm = numpy.bmat('New Wm')
        #~ Am_temp = Am*As[0]
        #~ Pb_temp = Wm*Wm.T + Am_temp*P*Am_temp.T
        #~ saveData(Eve.Obs,injectionTime[i+1],mbs[i],Pb_temp)  # Saves error bars ONLY ObsNum = 0
    #~ Am = Am*As[0]
    #~ # print 'Wm', Wm
    #~ # print 'Am', Am
    #~ # print 'P', P
    #~ Pb = Wm*Wm.T + Am*P*Am.T
    #~ return (mb, Pb, t1)

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
