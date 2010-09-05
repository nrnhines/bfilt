from neuron import h
import fitglobals
fitglobals.debugoff()
h.load_file('mulfit.hoc')
import nrnbfilt

import sys
n = sys.modules['nrnbfilt'].__str__()
n = n[n.find('/'):n.rfind('/')]
if len(n) and n[-1] != '/':
    n += '/'
h.load_file(n+'eonerunmlf.hoc')
del n

