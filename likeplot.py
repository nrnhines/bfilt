import numpy
import fitEKF
import os
import fitglobals
import copy

def calc1(MS,R0=None):
    if R0 == None:
        R0 = numpy.arange(0.05,1.0,0.05)
    LLmax = float('-inf')
    Mfit = copy.deepcopy(MS.M)
    print('Turning debugging OFF')
    fitglobals.debugoff()
    f = open('data.txt','w')
    f.write('# A0 LL\n')
    for A0 in R0:
        # Mfit.P.A = numpy.matrix([[A0, 0], [0,A1]])
        MS.setParams(A0,0,Mfit)
        LL = fitEKF.ekf(MS.Data,Mfit)
        if LL>LLmax:
            A0max = A0
            LLmax = LL
        f.write('%s %s\n' % (A0, LL))
        print A0, LL
    f.close()
    print 'Max Like =', LLmax, '@ A0 =', A0max
    simParams = MS.getParams()
    LLsim = MS.loglike(simParams[0],simParams[1])
    print 'Sim Like =', LLsim, '@ A0 =', simParams[0]

def calc2(MS,R0=None,R1=None):
    if R0 == None:
        R0 = numpy.arange(0.05,1.0,0.05)
    if R1 == None:
        R1 = numpy.arange(0.05,1.1,0.05)

    LLmax = float('-inf')
    Mfit = copy.deepcopy(MS.M)
    print('Turning debugging OFF')
    fitglobals.debugoff()
    f = open('data.txt','w')
    f.write('# A0 A1 LL\n')
    for A0 in R0:
        for A1 in R1:
            # Mfit.P.A = numpy.matrix([[A0, 0], [0,A1]])
            MS.setParams(A0,A1,Mfit)
            LL = fitEKF.ekf(MS.Data,Mfit)
            if LL>LLmax:
                A0max = A0
                A1max = A1
                LLmax = LL
            f.write('%s %s %s\n' % (A0, A1, LL))
            print A0, A1, LL
        f.write('\n')
        print ""
    f.close()
    print 'Max Like =', LLmax, '@ (A0,A1) =', [A0max, A1max]
    simParams = MS.getParams()
    LLsim = MS.loglike(simParams[0],simParams[1])
    print 'Sim Like =', LLsim, '@ (A0, A1) =', simParams[0], simParams[1]

    