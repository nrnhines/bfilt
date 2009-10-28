#import jacflow
import math
#import scipy
#import scipy.linalg
import numpy

def ekf(data, model):
    m0 = model.Initial
    P0 = model.P.InitialCov
    k = 0
    smll = 0
    
    # Main Filtering Loop
    while k < len(data):
        # Evaluate derivatives for prediction
        Times = model.FitEvents[k][0]
        ObsNum = model.FitEvents[k][1]
        mb = model.Sys.flow(Times, m0)
        Am = model.Sys.Dstate(Times, m0)
        Wm = model.Sys.Dnoise(Times, m0)
        
        # Prediction
        Pb = Wm*Wm.T + Am*P0*Am.T
        
        # Evaluate derivatives for update
        hh = model.Obs.meas(Times, mb, ObsNum)
        H = model.Obs.Dstate(Times, mb, ObsNum)
        V = model.Obs.Dnoise(Times, mb, ObsNum)
        
        # Update
        v = data[k] - hh
        S = H*Pb*H.T + V*V.T
        K = Pb*H.T*S.I
        P0 = Pb - K*S*K.T
        m0 = mb + K*v
        
        # Likelihood
        #mll = math.log(scipy.linalg.det(S)) + v.T*S.I*v
        mll = math.log(S) + v.T*S.I*v
        smll += mll
        k += 1
    # Keeps running sum of Minus-Log-Likelihood
    # But returns Positive Log Likeliood (for now)
    return -smll
