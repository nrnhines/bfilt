import hmExperiment
import numpy
import string

def conSeed2HR():
    rn = numpy.arange(.001,12.0001,.05).tolist()
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[2]],tstops=[[20]])
    F.save4plot("conSeed2HR",dict(tau01=rn,tau12=rn))
    print F.find()

def findMLSeed2():
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[2]],tstops=[[20]])
    F.find()
    return F.ML

def z2cSeed2HR():
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[2]],tstops=[[20]])
    F.find()
    fname = "/Users/seancarver/Code/neuron-dev/bfilt/conSeed2HR"
    f = open(fname+"_z.txt","r")
    LN = f.readline()
    LNS = string.split(LN)
    f.close()
    g = open(fname+"_c.txt","w")
    for zstr in LNS:
        z = float(zstr)
        c = F.cValue(z)
        g.write(str(c)+' ')
    g.close()
    
def readlikefile(fnamep):
    f = open(fnamep,"r")
    LN = f.readline()
    LNS = string.split(LN)
    xs = []
    for x in LNS:
        xs.append(float(x))
    f.close()
    return xs

def extendy(ymax=12.002):
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[2]],tstops=[[20]])
    F.find()
    fname = "/Users/seancarver/Code/neuron-dev/bfilt/conSeed2HR"
    xs = readlikefile(fname+"_x.txt")
    ys = readlikefile(fname+"_y.txt")
    zs = readlikefile(fname+"_z.txt")
    cs = readlikefile(fname+"_c.txt")
    dy = ys[-1] - ys[-2]
    f = open(fname+"NEW_y.txt", "w")
    for y in ys:
        f.write(str(y)+' ')
    y = ys[-1] + dy
    while y < ymax:
        f.write(str(y)+' ')
        y = y + dy
    f.close()
    f = open(fname+"NEW_x.txt", "w")
    for x in xs:
        f.write(str(x)+' ')
    i = 0
    p = {}
    p.update(F.guess)
    p.update(F.known)
    efun = hmExperiment.like4eval
    f = open(fname+"NEW_z.txt", "w")
    for x in xs:
        for y in ys:
            f.write(str(zs[i])+' ')
        y = ys[-1] + dy
        while y < ymax:
            p.update({"tau01":x,"tau12":y})
            newz = efun(p,F.structure,F.SysWData)
            f.write(str(newz)+' ')
             