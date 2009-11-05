from neuron import h
h.load_file('mulfit.hoc')
h.load_file('eonerunmlf.hoc')

if __name__ == '__main__':
  h.load_file('nrngui.hoc')
  from math import pi
  s = h.Section()
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
