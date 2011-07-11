import hmExperiment
import numpy

def conSeed2HR():
    rn = numpy.arange(.001,12.0001,.05).tolist()
    F = hmExperiment.fit(hmExperiment.ch3up,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3up,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[2]],tstops=[[20]])
    F.save4plot("conSeed2HR",dict(tau01=rn,tau12=rn))
    print F.find()
    
    