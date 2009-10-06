import jacflow
import math
import scipy
import scipy.linalg
import numpy
from neuron import h

R = numpy.matrix(1)

def ekf(data):
  k = 0
  smll = 0
  h.stdinit()
  P0 = numpy.matrix(numpy.eye(4))
 
  # Main Filtering Loop
  while k < len(data):

    # Evaluate derivatives for prediction
    L = jacflow.flowJac(data[k][0])
    mb = L[0]
    Am = L[1]
    Wm = L[2]
  
    # Prediction
    # Q identity; Pb = Wm*Q*Wm.T + Am*P0*Am.T
    Pb = Wm*Wm.T + Am*P0*Am.T

    # Evaluate derivatives for update
    L = jacflow.measJac()
    hh = L[0]
    H = L[1]
    V = L[2]

    # Update
    v = data[k][1] - hh
    # R identity; S = H*Pb*H.T + V*R*V.T
    print H.shape
    print Pb.shape
    print V.shape
    S = H*Pb*H.T + V*V.T
    K = Pb*H.T*S.I
    P0 = Pb - K*S*K.T
    m0 = mb + K*v
    h.cvode.yscatter(h.Vector(m0))
    h.cvode.re_init()

    # Marginal Measurement Likelihood
    smll += math.log(scipy.linalg.det(S)) + v.T*S.I*v

    k += 1
  return smll
