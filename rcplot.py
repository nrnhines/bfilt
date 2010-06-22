import pickle
import math
import numpy
import matplotlib
import matplotlib.pylab as pylab

def plot():
    (n,P,eP,A,eA) = numbers()

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
    (A,eA) = accuracy([T2,T8,T32,T128])
    return (n,P,eP,A,eA)

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
    return (A,eA) 
