from neuron import h
from math import pi
h.load_file('nrngui.hoc')
ss = [h.Section(), h.Section()]
for s in ss:
    s.diam = 10
    s.L = 100/pi/s.diam
    s.insert('pas')
    s.e_pas = 0
ss[0].g_pas = .5*.001
ss[1].g_pas = .1*.001
h.tstop = 20
h.dt = 0.1
h.v_init = 1
