import hmEnsemble

def fig(seed=0,x=[2.,4.],y=[2.,4.],nc0=5,nc1=5):
    M0 = hmEnsemble.ch3Ensemble(tau01=2,tau12=4,nchannels=nc0)
    M0.sim(seed=seed)
    L = []
    G = []
    for tau01 in x:
        L.append([])
        G.append([])
        for tau12 in y:
            M1 = hmEnsemble.ch3Ensemble(tau01=tau01,tau12=tau12,nchannels=nc1)
            L[-1].append(M1.likelihood(M0.simData))
            G[-1].append((tau01,tau12))
    return (L,G)
