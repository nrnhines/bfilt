import math
import numpy
import pylab
from neuron import h
import stochsim

sqrtEps = math.sqrt(numpy.finfo(numpy.double).eps)
sigp = stochsim.sigp
sigm = stochsim.sigm 
x = h.Vector()
value = h.Vector()
df = h.Vector()

def flowJac(noiseTimes):

  h.cvode.states(x)
  t0 = h.t
  print 't0' 
  print t0
  ss = [None]*(len(noiseTimes)+1)
  k = 0
  ss[0] = h.SaveState()
  ss[0].save()
  h.cvode.re_init()
  while k<len(noiseTimes):
    print noiseTimes[k]
    h.continuerun(noiseTimes[k])
    ss[k+1] = h.SaveState()
    ss[k+1].save()
    h.cvode.re_init()
    k += 1
  h.cvode.states(value)

  # Derivative with respect to state
  DFx = pylab.zeros((len(value),len(x)))
  k = 0
  while k<len(x):
    ss[0].restore(1)
    h.cvode.re_init()
    h.cvode.states(x)
    temp = x[k]
    if abs(temp) > 1:
      dx = sqrtEps*abs(temp)
    else:
      dx = sqrtEps
    x.x[k] = temp + dx
    h.cvode.yscatter(x)
    h.cvode.re_init()
    dx = x[k] - temp  # trick to reduce finite precision error
    h.continuerun(noiseTimes[len(noiseTimes)-1])
    h.cvode.states(df)
    df.sub(value)
    df.div(dx)
    DFx[:,k] = numpy.array(df);
    k += 1

  # Derivative with respect to noise
  DFn = pylab.zeros((len(value),len(noiseTimes)))
  noiseTimes.insert(0,t0)
  k = 0
  while k<len(noiseTimes)-1:
    ss[k+1].restore(1)
    h.cvode.re_init()
    if k < len(noiseTimes)-2:
      dx = sqrtEps
      h.cvode.states(x)
      x.x[0] += sigp*math.sqrt(noiseTimes[k+1]-noiseTimes[k])*dx
      h.cvode.yscatter(x)
      h.cvode.re_init()
      h.continuerun(noiseTimes[len(noiseTimes)-1])
      h.cvode.states(df)
      df.sub(value)
      df.div(dx)
      DFn[:,k] = numpy.array(df);
    else:
      DFn[:,k] = numpy.array([sigp*math.sqrt(noiseTimes[k+1]-noiseTimes[k]),0,0,0])
    k += 1

  return [numpy.matrix(value).T, numpy.matrix(DFx), numpy.matrix(DFn)]

def measJac():
  x = h.Vector();
  h.cvode.states(x)
  return [numpy.matrix(x.x[0]),numpy.matrix([1,0,0,0]),numpy.matrix(sigm)]

def flow1Test(a, x, dt, sp):
  return  [numpy.matrix(x*math.exp(-a*dt)), numpy.matrix(math.exp(-a*dt)), numpy.matrix(sp*math.sqrt(dt))]
  
def measTest(x, sm):
  return [numpy.matrix(x), numpy.matrix(1), numpy.matrix(sm)]
