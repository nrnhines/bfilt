import hmlike
import pickle

def getp(seed):
    print "seed", seed
    K = hmlike.HML()
    K.sim(seed)
    return (K.p_true(),K.gopt,K.mle)

def runreps(reps=5,ntraj=16):
    print "reps",reps
    print "ntraj",ntraj
    ss = []
    ps = []
    gs = []
    xs = []
    multiplier = 100000
    for i in range(reps):
        seeds = []
        for trajnum in range(ntraj):
            seeds.append(i+trajnum*multiplier)
        (p,g,x) = getp(seeds)
        ss.append(seeds)
        ps.append(p)
        gs.append(g)
        xs.append(x)
        print 'gopt', g
    fname = "reps"+str(reps)+"ntraj"+str(ntraj)+".pkl"
    f = open(fname,'w')
    pickle.dump(ss,f)
    pickle.dump(ps,f)
    pickle.dump(gs,f)
    pickle.dump(xs,f)
    f.close
