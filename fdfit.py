import jacflow
import math
from myscypi import linalg
import numpy
from neuron import h

R = numpy.matrix(1)

def ekf(data, aparam):
  sigmap = 0.01
  sigmam = 0.001
  dt = 0.01
  m0 = numpy.matrix(1)
  P0 = numpy.matrix(1)

  k = 0
  smll = 0
  # h.stdinit()
  # P0 = numpy.matrix(numpy.eye(4))
 
  # Main Filtering Loop
  while k < len(data):

    # Evaluate derivatives for prediction
    # L = jacflow.flowJac(data[k][0])
    L = jacflow.flow1Test(aparam, m0, dt, sigmap)
    mb = L[0]
    Am = L[1]
    Wm = L[2]
  
    # Prediction
    # Q identity; Pb = Wm*Q*Wm.T + Am*P0*Am.T
    Pb = Wm*Wm.T + Am*P0*Am.T
    # print Wm.shape
    # print Am.shape
    # print('Pb')
    # print linalg.eig(Pb) 

    # Evaluate derivatives for update
    # L = jacflow.measJac()
    L = jacflow.measTest(mb, sigmam)
    hh = L[0]
    H = L[1]
    V = L[2]

    # Update
    v = data[k][1] - hh
    # print 'v'
    # print v
    # R identity; S = H*Pb*H.T + V*R*V.T
    # print H.shape
    # print Pb.shape
    # print V.shape
    # print 'S' 
    S = H*Pb*H.T + V*V.T
    # print S 
    K = Pb*H.T*S.I
    P0 = Pb - K*S*K.T
    # print 'P0' 
    # print linalg.eig(P0) 
    m0 = mb + K*v
    # h.cvode.yscatter(h.Vector(m0))
    # h.cvode.re_init()

    # Marginal Measurement Likelihood
    mll = math.log(linalg.det(S)) + v.T*S.I*v
    # print 'mll'
    # print mll
    smll += mll
    
    # print linalg.eig(S) 
    k += 1
    spll = -smll
  return spll
