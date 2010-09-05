from neuron import h
cv = h.cvode

use_fixed_step = False
class FSPanel(object):
  def __init__(self):
    self.use_fixed_step = 0.0
  def panel(self):
    h.xpanel("Variable or Fixed step")
    h.xcheckbox('use_fixed_step', (self,'use_fixed_step'))
    h.xpanel()

fs = FSPanel()
#fs.panel() # seems to break testcr.onerun() when put here.

def active(x):
  return cv.active(x)

def  atol(x):
  return cv.atol(x)

def f(t, v, vp):
  return cv.f(t, v, vp)

def re_init():
  return cv.re_init()

def solve(x):
  global fs
  #print 'before solve %g t=%g'%(x,h.t)
  if (fs.use_fixed_step == 1.0):
    h.dt = .025
    while (h.t < (x - h.float_epsilon)):
      cv.fixed_step()
    re_init()
  else:
    cv.solve(x)
  #print 'after solve %g t=%g'%(x,h.t)
  return 1.

def statename(i, sref, j):
  return cv.statename(i, sref, j)

def states(v):
  return cv.states(v)

def yscatter(v):
  cv.yscatter(v)
  return 1.0

