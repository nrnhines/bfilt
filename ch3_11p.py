from neuron import h
import fitglobals
fitglobals.debugoff()
h.load_file('mulfit.hoc')
h.load_file('eonerunmlf.hoc')
import nrnbfilt
h.load_file('ch3_11p.ses')

h('objref nb')
h.nb = h.List("PythonObject").o(0)
