from neuron import h
from scipy.optimize import fmin_bfgs
import numpy

h('''
mulfit_optimizers_append("BFGS", "BFGSWrap")
begintemplate BFGSWrap
public prun, showopt, save_optimizer, bfgs
objref bfgs
proc init() { bfgs=new PythonObject()  bfgs=bfgs.bfgswrap.BFGS($o1) }
func prun() { return bfgs.prun($o1) }
proc showopt() {}
proc save_optimizer() {}
endtemplate BFGSWrap
'''
)

class BFGS(object):

  def __init__(self, pf):
    print 'created BFGS instance'
    self.pf = pf
    self.pf.def_parmlist_use()
    self.opt = None
    self.parmvec = h.Vector()

  def prun(self, opt):
    self.opt = opt
    x0 = numpy.array(opt.start)
    print x0
    e = fmin_bfgs(self.efun, x0, epsilon=1e-4)
    print e
    return e[1]

  #callback for fmin_bfgs
  def efun(self, x):
    v = self.parmvec.from_python(x)
    if self.opt:
      e = self.opt.optwrap_efun(v.size(), v._ref_x[0])
    else:
      e = self.pf.efun(v.size(), v._ref_x[0])
    return e
