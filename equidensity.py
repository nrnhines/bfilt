import svd
import EKF
import math
import numpy

def equipoints(k,state0,state1,num):
    mat = [[EKF.Ps[k][state0, state0],EKF.Ps[k][state0, state1]],[EKF.Ps[k][state1, state0],EKF.Ps[k][state1, state1]]]
    decomp = svd.svd(mat)
    U = decomp[0]
    std0 = math.sqrt(decomp[1][0])
    std1 = math.sqrt(decomp[1][1])
    SLam = [[std0,0.0],[0.0,std1]]
    m = [[EKF.ms[k][state0,0]],[EKF.ms[k][state1,0]]]
    m = numpy.matrix(m)
    SLam = numpy.matrix(SLam)
    U = numpy.matrix(U)
    x = []
    y = []
    for angle in numpy.arange(0,2*math.pi+2*math.pi/num,2*math.pi/num):
        circlepoint = numpy.matrix([[math.cos(angle)],[math.sin(angle)]])
        xypoints = m + U*SLam*circlepoint
        x.append(xypoints[0,0])
        y.append(xypoints[1,0])
    return (x,y)
