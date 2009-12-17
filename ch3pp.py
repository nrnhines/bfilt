from neuron import h
import fitglobals
fitglobals.debugoff()
h.load_file('mulfit.hoc')
h.load_file('eonerunmlf.hoc')
import nrnbfilt
h.load_file('ch3pp.ses')
