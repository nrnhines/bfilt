from neuron import h

class Flow(object):
  def __init__(self):
    self.ss = h.SaveState()
    self.ss.save()
    self.t0 = h.t

  def flow(self, x0, t1, x1):
    self.ss.restore(1)
    h.cvode.yscatter(x0)
    h.cvode.re_init()
    h.initPlot()
    h.continuerun(t1) # or h.cvode.solve(t1)
    h.cvode.states(x1)
    return x1

if __name__ == '__main__':
  h.load_file('hh.ses')
  h.Graph[0].exec_menu('Keep Lines')
  h.stdinit()
  h.continuerun(1.5)
  flow = Flow()
  x = h.Vector()
  x1 = h.Vector()
  h.cvode.states(x)
  for i in range(5):
    print i, x
    x.x[0] += 5
    flow.flow(x, h.tstop, x1)
