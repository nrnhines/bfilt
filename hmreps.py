import hmlike
import pickle

def getp(seed):
    print "seed", seed
    K = hmlike.HML()
    K.sim(seed)
    return K.p_true()

def reps(n):
    ps = []
    for seed in range(n):
        p = getp(seed)
        ps.append(p)
    fname = "reps"+str(n)+".pkl"
    f = open(fname,'w')
    pickle.dump(ps,f)
    f.close
