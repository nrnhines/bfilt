import hmlike
import pickle

def getp(seed):
    print "seed", seed
    K = hmlike.HML()
    K.sim(seed)
    return (K.p_true(),K.gopt,K.mle)

def reps(n):
    ps = []
    gs = []
    for seed in range(n):
        (p,g,x) = getp(seed)
        ps.append(p)
        gs.append(g)
        xs.append(x)
        print 'gtol', g
    fname = "reps"+str(n)+".pkl"
    f = open(fname,'w')
    pickle.dump(ps,f)
    pickle.dump(gs,f)
    f.close
