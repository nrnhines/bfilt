import pickle
import math
import numpy
import matplotlib
import pylab

def plot():
    (n,P,eP,A,eA,ALines) = numbers()
    pylab.subplot(2,1,1)
    pylab.errorbar(n,P,eP,ecolor='r',linewidth=3,fmt=None)
    pylab.semilogx()
    pylab.xticks([2,8,32,128],['2','8','32','128'])
    pylab.axis([1,256,4,18])
    pylab.ylabel('Precision of 95% Confidence')
    pylab.title('Estimating Time Constant of RC circuit, Gaussian Noise, 100 trials')
    pylab.subplot(2,1,2)
    pylab.plot(n,A,'*',color='r',markersize=12)
    #pylab.axhspan(ALines[0],ALines[0],color='b')
    pylab.axhspan(ALines[1],ALines[1],color='b')
    pylab.axhspan(ALines[2],ALines[2],color='b')
    pylab.semilogx()
    pylab.xticks([2,8,32,128],['2','8','32','128'])
    pylab.axis([1,256,0.9,1.0])
    pylab.text(1.5,0.95,'Accurate')
    pylab.text(1.5,0.99,'Under-Confident')
    pylab.text(1.5,0.91,'Over-Confident')
    pylab.xlabel('Number of Data Points')
    pylab.ylabel('Accuracy of 95% Confidence')
    pylab.ion()

def numbers():
    f = open("T2pickle.pkl","r")
    T2 = pickle.load(f)
    f.close()
    f = open("T8pickle.pkl","r")
    T8 = pickle.load(f)
    f.close()
    f = open("T32pickle.pkl","r")
    T32 = pickle.load(f)
    f.close()
    f = open("T128pickle.pkl","r")
    T128 = pickle.load(f)
    f.close()
    n = [2,8,32,128]
    (P,eP) = precision([T2,T8,T32,T128])
    (A,eA,ALines) = accuracy([T2,T8,T32,T128])
    return (n,P,eP,A,eA,ALines)

def precision(L):
    P = []
    eP = []
    for l in L:
        ps = []
        ntot = 0
        for t in l:
            ntot += 1
            ps.append(t.precision)
        P.append(numpy.mean(ps))
        eP.append(numpy.std(ps)/math.sqrt(float(ntot)))
    return (P,eP)


def accuracy(L):
    A = []
    eA = []
    for l in L:
        nin = 0
        ntot = 0
        for t in l:
            ntot += 1
            if t.covers:
                nin += 1
        p = (float(nin)/float(ntot))
        A.append(p)
        eA.append(math.sqrt(p*(1-p)/float(ntot)))
        SE95 = math.sqrt(0.95*0.05/float(ntot))
        ALines = [0.5,0.95-SE95,0.95+SE95]
    return (A,eA,ALines)
