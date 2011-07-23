import hmExperiment
import numpy
import string

# Keep this function, used for plots
def x5c1(seeds=[[12,22,32,42,52]],fname="conX5c1SX2cor1cF"):
    rn = numpy.arange(.001,16.51,.05).tolist()
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=1))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=1),seeds=seeds,tstops=[[20]])
    F.save4plot(fname,dict(tau01=rn,tau12=rn,nchannels=1))
    print F.find()

def x1c1():
    x5c1(seeds=[[2]],fname="conX1c1S2cor1cF")

def conV5(Vchar=5,fname="conV5",seeds=[[2]]):
    rn = numpy.arange(.001,16.01,.05).tolist()
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(Vchar01=Vchar,Vchar12=Vchar,nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5,Vchar01=Vchar,Vchar12=Vchar),seeds=seeds,tstops=[[20]])
    F.save4plot(fname,dict(tau01=rn,tau12=rn,nchannels=5,Vchar01=Vchar,Vchar12=Vchar))
    print F.find()

def conV10():
    Vchar=10
    fname="conV10"
    conV5(Vchar=Vchar,fname=fname);

def conV20(seeds=[[2]]):
    Vchar=20
    fname="conV20"
    conV5(Vchar=Vchar,fname=fname,seeds=seeds);

def conV40():
    Vchar=40
    fname="conV40"
    conV5(Vchar=Vchar,fname=fname);

def conSHR16(seeds=[[2]],fname="conSeed2HR16"):
    rn = numpy.arange(.001,16.01,.05).tolist()
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=seeds,tstops=[[20]])
    F.save4plot(fname,dict(tau01=rn,tau12=rn,nchannels=5))
    print F.find()

def conSeed0HR16():
    conSHR16([[0]],"conSeed0HR16")

def conSeed1HR16():
    conSHR16([[1]],"conSeed1HR16")

def conSeed2HR16():
    conSHR16([[2]],"conSeed2HR16")

def conSeed3HR16():
    conSHR16([[3]],"conSeed3HR16")

def conSeed4HR16():
    conSHR16([[4]],"conSeed4HR16")

def conSeed5HR16():
    conSHR16([[5]],"conSeed5HR16")

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

def extendy(ymax=13.002):
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[2]],tstops=[[20]])
    F.find()
    fname = "conSeed2HR"
    xs = readlikefile(fname+"_x.txt")
    ys = readlikefile(fname+"_y.txt")
    zs = readlikefile(fname+"_z.txt")
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
            print "tau01", x, "tau12", y
            newz = efun(p,F.structure,F.SysWData)
            f.write(str(newz)+' ')
            y = y + dy 
