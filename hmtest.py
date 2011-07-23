import hmm
import hmExperiment
import numpy
import hmRun

def chaintest():
    HC1 = hmm.ch3chain([-65.,20.,-65.])
    HC2 = hmm.ch3chain([-65.,20.,-65.],tau01=10,tau12=15)
    HC1.sim(seeds=[1,2,3,4,5],tstops=[20.,20.])
    HC2.sim(seeds=[6,7,8,9,10],tstops=[20.,20.])
    l11 = HC1.likelihood(HC1)
    l12 = HC1.likelihood(HC2)
    l21 = HC2.likelihood(HC1)
    l22 = HC2.likelihood(HC2)
    assert(False)

def experchntest():
    rn = numpy.arange(2.,5.,2.).tolist()
    F = hmExperiment.fit(hmExperiment.ch3chn,dict(tau01=2.,tau12=4.),dict(nchannels=5))
    F.sim(hmExperiment.ch3chn,dict(tau01=2.,tau12=4.,nchannels=5),seeds=[[1,2]],dt=.1,tstops=[[20,20]])
    F.save4plot("testLR",dict(tau01=rn,tau12=rn))
 
def V20test():    
    hmRun.conV20()

V20test()
#chaintest()
#experchntest()
