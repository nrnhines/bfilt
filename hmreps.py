import hmlike
import pickle

def getp(seed):
    print "seed", seed
    K = hmlike.HML()
    K.sim(seed)
    return (K.p_true(),K.gtol)

def reps(n):
    ps = []
    gs = []
    for seed in range(n):
        (p,g) = getp(seed)
        ps.append(p)
        gs.append(g)
        print 'gtol', g
    fname = "reps"+str(n)+".pkl"
    f = open(fname,'w')
    pickle.dump(ps,f)
    pickle.dump(gs,f)
    f.close
