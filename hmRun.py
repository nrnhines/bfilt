import hmExperiment
import numpy

def conSeed2():
    rn = numpy.arange(.01,15.001,.5).tolist()
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[2]],tstops=[[20]])
    F.save4plot("conSeed2",dict(tau01=rn,tau12=rn))
    print F.find()
    
    