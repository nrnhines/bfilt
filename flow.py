from neuron import h
import cvodewrap

class Flow(object):
    def __init__(self):
        self.ss = h.SaveState()
        self.ss.save()
        self.t0 = h.t
    
    def flow(self, x0, t1, x1):
        self.ss.restore(1)
        cvodewrap.yscatter(x0)
        cvodewrap.re_init()
        h.initPlot()
        h.continuerun(t1) # or cvodewrap.solve(t1)
        cvodewrap.states(x1)
        return x1

if __name__ == '__main__':
    h.load_file('hh.ses')
    h.Graph[0].exec_menu('Keep Lines')
    h.stdinit()
    h.continuerun(1.5)
    flow = Flow()
    x = h.Vector()
    x1 = h.Vector()
    cvodewrap.states(x)
    for i in range(5):
        print i, x
        x.x[0] += 5
        flow.flow(x, h.tstop, x1)
