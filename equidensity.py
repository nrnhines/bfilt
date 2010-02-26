import svd
import EKF
import math
import numpy
import Gnuplot
import nrnbfilt

G = Gnuplot.Gnuplot()

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
    xy = []
    for angle in numpy.arange(0,2*math.pi+2*math.pi/num,2*math.pi/num):
        circlepoint = numpy.matrix([[math.cos(angle)],[math.sin(angle)]])
        xypoints = m + U*SLam*circlepoint
        xy.append([xypoints[0,0], xypoints[1,0]])
    return xy

def frame(k,state0,state1,num):
    xy = equipoints(k,state0,state1,num)
    global G
    G.plot(xy)
    G.reset()

def movie(state0,state1,num=200):
    n = len(EKF.ms)
    for k in range(n):
        print k
        frame(k,state0,state1,num)
        raw_input()