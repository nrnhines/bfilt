import hmm

chaintest()

def chaintest():
    HC = hmm.ch3chain([-65.,20.,-65.])
    HC.sim(seeds=[1,2,3,4,5],tstops=[20.,20.])
    HC.likelihood(HC)
    print "Successful..."
    