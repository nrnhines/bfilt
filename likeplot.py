import runTwoDecay
import numpy
import fitEKF

A0range = numpy.arange(0.05,1.0,0.05)
A1range = numpy.arange(0.05,1.1,0.05)

f = open('data.txt','w')
f.write('# A0 A1 LL\n')
for A0 in A0range:
    for A1 in A1range:
        runTwoDecay.M.P.A = numpy.matrix([[A0, 0], [0,A1]])
        LL = fitEKF.ekf(runTwoDecay.Data,runTwoDecay.M)
        f.write('%s %s %s\n' % (A0, A1, LL.tolist()[0][0]))
        print A0, A1, LL.tolist()[0][0]
    f.write('\n')
    print ""
