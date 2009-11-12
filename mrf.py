from neuron import h
import fitglobals
h.load_file('mulfit.hoc')
h.load_file('eonerunmlf.hoc')
import nrnbfilt

def test():
  global s, h
  fitglobals.debugon()
  h.load_file('nrngui.hoc')
  from math import pi
  h('create soma')
  s = h.soma
  s.diam = 10
  s.L = 100/pi/s.diam
  s.insert('pas')
  s.e_pas = 0
  s.g_pas = .5*.001
  h.tstop = 10
  h.dt = 0.1
  h.v_init = 1
  h.load_file('pasfit.ses')
  h.cvode_active(1)
  #import testext # if launching python instead of nrniv

if __name__ == '__main__':
  test()
